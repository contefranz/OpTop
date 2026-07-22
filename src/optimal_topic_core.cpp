#include <RcppArmadillo.h>
// [Rcpp::depends(RcppArmadillo)]

#include <algorithm>
#include <numeric>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Model-agnostic core of the Test 1 statistic (Eq. 8 in Lewis & Grossetti
// 2022): the caller adapts whatever topic model it holds to a document-topic
// matrix theta (J_model x K, rows summing to 1) and a topic-word matrix phi
// (K x W, rows summing to 1) and calls this once per model. The hot path is
// free of R API objects (no S4 slot access), which is what allows the
// per-document loop to run under OpenMP.
//
// Threading contract: documents are independent, so within each block the
// per-document envelope/Pearson work runs in parallel, each document writing
// its own result slot; the slots are then accumulated serially in document
// order. The floating-point summation order (and the exported envelope
// layout) therefore never depends on the schedule, and the output is
// bit-identical for any n_threads. The block densification of the sparse
// dfm runs under the same contract (one column per document).
//
// Tie handling: envelope order is defined by the comparator
// (probability descending, original index ascending), the semantics of R's
// stable order(X, decreasing = TRUE). Ties in the fitted probabilities
// (common for count-based estimators such as collapsed Gibbs) therefore
// resolve identically on every platform and match the R reference
// implementation exactly.
//
// Envelope selection: the statistic needs the smallest descending-sorted
// head whose cumulative mass exceeds q, so a full O(W log W) sort is not
// required. The head is located with std::nth_element on a candidate
// boundary that doubles until the head mass crosses q; each widening sorts
// only the new slice (the partition property guarantees slices are already
// ordered relative to each other), for an overall cost of
// O(W · widenings + P_j log P_j) instead of O(W log W). The initial
// boundary is a fixed 1024 (clamped to W): starting at a fraction of W
// would sort thousands of entries per document at large vocabularies when
// the typical P_j is tens to hundreds, and because the comparator is a
// strict total order the boundary trajectory cannot change the selected
// envelope, only the work done to find it.

namespace {

struct EnvelopeCmp {
    const double* x;
    bool operator()(const arma::uword a, const arma::uword b) const {
        if (x[a] != x[b]) return x[a] > x[b];
        return a < b;
    }
};

}  // namespace

//' @keywords internal
// [[Rcpp::export]]
Rcpp::List optimal_topic_core(const arma::mat& theta,
                              const arma::mat& phi,
                              const Rcpp::IntegerVector& dfm_p,
                              const Rcpp::IntegerVector& dfm_i,
                              const Rcpp::NumericVector& dfm_x,
                              int n_terms_in,
                              double q,
                              const arma::uvec& doc_map,
                              bool return_envelope,
                              int n_threads)
{
    // The weighted dfm arrives transposed once by the caller (W x J) as the
    // raw CSC slots of the dgCMatrix (dfm_p column pointers, dfm_i row
    // indices, dfm_x values): a document is a contiguous CSC column, and the
    // Rcpp vectors are zero-copy views of R memory, so no Armadillo sparse
    // container is ever constructed. This removes both the per-call copy of
    // the corpus and the J x W addressability check of SpMat::init(), so the
    // shape of the corpus is unbounded; only the nonzero count is bounded
    // (by R's dgCMatrix container, at 2^31 - 1).
    if (dfm_p.size() < 1) {
        Rcpp::stop("the dfm column pointers are empty");
    }
    if (n_terms_in < 0) {
        Rcpp::stop("the vocabulary size must be nonnegative");
    }
    const arma::uword n_docs = static_cast<arma::uword>(dfm_p.size() - 1);
    const arma::uword n_terms = static_cast<arma::uword>(n_terms_in);
    const arma::uword current_k = theta.n_cols;

    if (dfm_i.size() != dfm_x.size()) {
        Rcpp::stop("the dfm index and value slots disagree in length");
    }
    if (static_cast<R_xlen_t>(dfm_p[dfm_p.size() - 1]) != dfm_i.size()) {
        Rcpp::stop("the dfm column pointers do not close on the nonzeros");
    }
    if (phi.n_rows != current_k) {
        Rcpp::stop("theta and phi disagree on the number of topics");
    }
    if (phi.n_cols != n_terms) {
        Rcpp::stop("phi and the dfm disagree on the vocabulary size");
    }
    if (doc_map.n_elem != n_docs) {
        Rcpp::stop("the document mapping must have one entry per document");
    }
    if (n_docs == 0) {
        Rcpp::stop("the dfm has no documents");
    }
    if (doc_map.max() >= theta.n_rows) {
        Rcpp::stop("document mapping points past the rows of theta: "
                   "the dfm and the model are misaligned");
    }
    if (n_threads < 1) {
        n_threads = 1;
    }

    // raw pointers hoisted once: safe to read inside OpenMP regions
    const int* dfm_pp = dfm_p.begin();
    const int* dfm_ii = dfm_i.begin();
    const double* dfm_xx = dfm_x.begin();

    // documents are processed in fixed-size blocks: one BLAS product per
    // block replaces a dense W x K temporary per document, while keeping
    // memory bounded for large corpora
    const arma::uword block_size = 256;

    // W x K copy so each block's fitted probabilities are one gemm
    const arma::mat tww = phi.t();

    double sum_chi = 0.0;
    double sum_df = 0.0;

    // per-document envelope, exported for the calibration layer: the bin
    // structure (kept-word probabilities plus the collapsed min-bin mass)
    // depends only on the model, so the null distribution of the statistic
    // can be simulated or moment-matched from these bins alone
    std::vector<double> bin_probs;
    std::vector<int> bin_counts;
    if (return_envelope) {
        bin_probs.reserve(n_docs * 128);
        bin_counts.reserve(n_docs);
    }

    // per-document slots for one block, written in parallel and reduced
    // serially in document order
    std::vector<double> chi_slot(block_size);
    std::vector<double> df_slot(block_size);
    std::vector<std::vector<double>> env_slot(return_envelope ? block_size : 0);

    // block buffers allocated once: at large W the fitted-probability block
    // is on the order of 100 MB, so a fresh allocation per block would put
    // the allocator on the hot path. The ragged final block writes and
    // reads only its first block_len columns; stale columns are never
    // touched.
    arma::mat dfm_block(n_terms, block_size);
    arma::mat X_block(n_terms, block_size);
    arma::mat theta_t_blk(current_k, block_size);

    for (arma::uword block_start = 0; block_start < n_docs; block_start += block_size)
    {
        Rcpp::checkUserInterrupt();

        const arma::uword block_end = std::min(block_start + block_size, n_docs) - 1;
        const int block_len = static_cast<int>(block_end - block_start + 1);

        // gather the block's theta rows transposed (K x block_len), then one
        // gemm straight into the preallocated block through no-copy views:
        // no per-block gather, transpose, or product temporaries
        for (int j_col = 0; j_col < block_len; ++j_col) {
            const arma::uword r = doc_map[block_start + static_cast<arma::uword>(j_col)];
            double* dst = theta_t_blk.colptr(j_col);
            for (arma::uword k = 0; k < current_k; ++k) {
                dst[k] = theta.at(r, k);
            }
        }
        {
            const arma::uword bl = static_cast<arma::uword>(block_len);
            arma::mat X_view(X_block.memptr(), n_terms, bl, false, true);
            const arma::mat T_view(theta_t_blk.memptr(), current_k, bl,
                                   false, true);
            // no operand aliases the destination, so Armadillo evaluates
            // the product directly into the external memory
            X_view = tww * T_view;
        }

        // densify the block, one sparse column per document, in parallel
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
        for (int j_col = 0; j_col < block_len; ++j_col) {
            double* col = dfm_block.colptr(j_col);
            std::fill(col, col + n_terms, 0.0);
            const arma::uword j_doc = block_start + static_cast<arma::uword>(j_col);
            const int c0 = dfm_pp[j_doc];
            const int c1 = dfm_pp[j_doc + 1];
            for (int idx = c0; idx < c1; ++idx) {
                col[dfm_ii[idx]] = dfm_xx[idx];
            }
        }

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
        {
            // per-thread work buffers, reused across the documents a thread
            // processes
            std::vector<arma::uword> idx(n_terms);
            std::vector<double> x_buf(n_terms);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (int j_col = 0; j_col < block_len; ++j_col)
            {
                const double* x_src = X_block.colptr(j_col);
                std::copy(x_src, x_src + n_terms, x_buf.begin());
                const double* x = x_buf.data();
                const double* o = dfm_block.colptr(j_col);

                std::iota(idx.begin(), idx.end(), arma::uword(0));
                const EnvelopeCmp cmp{x};

                // locate the smallest sorted head whose cumulative mass
                // strictly exceeds q = 1 - I^K (Eq. 8; the crossing word is
                // kept, so the collapsed tail keeps mass strictly below I^K,
                // cf. footnote 5 of the paper): widen a selection boundary
                // until the head mass crosses q, sorting only new slices
                arma::uword sorted_upto = 0;
                arma::uword boundary = std::min<arma::uword>(n_terms, 1024);
                double cum = 0.0;
                arma::uword p_j = 0;
                bool crossed = false;

                while (!crossed) {
                    if (boundary > sorted_upto) {
                        if (boundary < n_terms) {
                            std::nth_element(idx.begin() + sorted_upto,
                                             idx.begin() + boundary,
                                             idx.end(), cmp);
                        }
                        std::sort(idx.begin() + sorted_upto,
                                  idx.begin() + boundary, cmp);
                    }
                    while (sorted_upto < boundary) {
                        cum += x[idx[sorted_upto]];
                        ++sorted_upto;
                        if (cum > q) {
                            p_j = sorted_upto;
                            crossed = true;
                            break;
                        }
                    }
                    if (!crossed) {
                        if (boundary >= n_terms) {
                            // the mass never crosses q: every word is kept
                            p_j = n_terms;
                            break;
                        }
                        boundary = std::min<arma::uword>(n_terms, boundary * 2);
                    }
                }

                const arma::uword n_tail = n_terms - p_j;

                // Pearson over the head, accumulated in envelope order
                double pearson = 0.0;
                double head_x = 0.0;
                double head_o = 0.0;
                for (arma::uword p = 0; p < p_j; ++p) {
                    const double xp = x[idx[p]];
                    const double op = o[idx[p]];
                    const double diff = op - xp;
                    pearson += diff * diff / xp;
                    head_x += xp;
                    head_o += op;
                }

                std::size_t n_bins = p_j;
                double tail_x = 0.0;
                if (n_tail > 0) {
                    // tail sums via the document totals: O(W) once instead
                    // of a second pass over the unsorted remainder
                    double total_x = 0.0;
                    double total_o = 0.0;
                    for (arma::uword w = 0; w < n_terms; ++w) {
                        total_x += x[w];
                        total_o += o[w];
                    }
                    tail_x = total_x - head_x;
                    const double tail_diff = (total_o - head_o) - tail_x;
                    pearson += tail_diff * tail_diff / tail_x;
                    n_bins += 1;
                }

                chi_slot[j_col] = double(n_bins) * pearson;
                df_slot[j_col] = double(n_bins - 1);

                if (return_envelope) {
                    std::vector<double>& bins = env_slot[j_col];
                    bins.clear();
                    bins.reserve(n_bins);
                    for (arma::uword p = 0; p < p_j; ++p) {
                        bins.push_back(x[idx[p]]);
                    }
                    if (n_tail > 0) {
                        bins.push_back(tail_x);
                    }
                }
            }
        }

        // serial reduction in document order: the summation order (and the
        // envelope layout) is independent of the schedule
        for (int j_col = 0; j_col < block_len; ++j_col) {
            sum_chi += chi_slot[j_col];
            sum_df += df_slot[j_col];
            if (return_envelope) {
                const std::vector<double>& bins = env_slot[j_col];
                bin_probs.insert(bin_probs.end(), bins.begin(), bins.end());
                bin_counts.push_back(int(bins.size()));
            }
        }
    }

    // one row: the number of topics (the "topic" column consumed by the
    // caller), the raw Test 1 statistic and its degrees of freedom
    // sum_j P_j; standardization and p-values are computed by the R caller,
    // where the tail convention is visible and testable
    arma::mat Chi_K(1, 3);
    Chi_K.row(0) = arma::rowvec({
        double(current_k),
        sum_chi,
        sum_df
    });

    if (return_envelope) {
        return Rcpp::List::create(
            Rcpp::Named("stat") = Chi_K,
            Rcpp::Named("bin_probs") = Rcpp::NumericVector(bin_probs.begin(), bin_probs.end()),
            Rcpp::Named("bin_counts") = Rcpp::IntegerVector(bin_counts.begin(), bin_counts.end())
        );
    }
    return Rcpp::List::create(Rcpp::Named("stat") = Chi_K);
}
