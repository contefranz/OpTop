#include <RcppArmadillo.h>
// [Rcpp::depends(RcppArmadillo)]

// Model-agnostic core of the Test 1 statistic (Eq. 8 in Lewis & Grossetti
// 2022): the caller adapts whatever topic model it holds to a document-topic
// matrix theta (J_model x K, rows summing to 1) and a topic-word matrix phi
// (K x W, rows summing to 1) and calls this once per model. Keeping the hot
// path free of R API objects (no S4 slot access) is also what makes it
// eligible for OpenMP parallelization over documents later on.

//' @keywords internal
// [[Rcpp::export]]
Rcpp::List optimal_topic_core(const arma::mat& theta,
                              const arma::mat& phi,
                              const arma::sp_mat& dfm_t,
                              double q,
                              const arma::uvec& doc_map,
                              bool return_envelope)
{
    // dfm_t is the weighted dfm transposed once by the caller (W x J): a
    // document is a contiguous CSC column, so densifying a block of
    // documents is a cheap column-range extraction instead of
    // element-by-element sparse row reads
    const arma::uword n_docs = dfm_t.n_cols;
    const arma::uword n_terms = dfm_t.n_rows;
    const arma::uword current_k = theta.n_cols;

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

    for (arma::uword block_start = 0; block_start < n_docs; block_start += block_size)
    {
        Rcpp::checkUserInterrupt();

        const arma::uword block_end = std::min(block_start + block_size, n_docs) - 1;

        // W x b slabs, one column per document
        const arma::mat dfm_block(dfm_t.cols(block_start, block_end));
        const arma::uvec theta_rows = doc_map.subvec(block_start, block_end);
        const arma::mat X_block = tww * arma::mat(theta.rows(theta_rows)).t();

        for (arma::uword j_col = 0; j_col < X_block.n_cols; ++j_col)
        {
            // get cumulative sum of sorted word probabilities
            arma::vec X = X_block.col(j_col);
            arma::vec weighted_dfm_j_doc = dfm_block.col(j_col);

            arma::uvec sort_indices = arma::sort_index(X, "descend");
            X = X(sort_indices);
            arma::vec X_cumsum = arma::cumsum(X);
            weighted_dfm_j_doc = weighted_dfm_j_doc(sort_indices);

            // Eq. (8): the P_j important words are the smallest head whose
            // cumulative mass strictly exceeds q = 1 - I^K (the crossing
            // word is kept, so the collapsed tail keeps mass strictly below
            // I^K, cf. footnote 5), and the document contributes
            // (bins) * Pearson with df = bins - 1, where bins = P_j + 1
            // with the collapsed min bin
            auto first_greater_than_q = std::find_if(X_cumsum.begin(),
                                                     X_cumsum.end(),
                                                     [&q](double val) {return val > q;});
            const std::size_t p_j = (first_greater_than_q == X_cumsum.end())
                ? X.n_elem
                : std::distance(X_cumsum.begin(), first_greater_than_q) + 1;
            const std::size_t n_tail = X.n_elem - p_j;

            const arma::vec head_diff = weighted_dfm_j_doc.head(p_j) - X.head(p_j);
            double pearson = arma::sum(head_diff % head_diff / X.head(p_j));

            std::size_t n_bins = p_j;
            double tail_X = 0.0;
            if (n_tail > 0) {
                tail_X = arma::sum(X.tail(n_tail));
                const double tail_diff = arma::sum(weighted_dfm_j_doc.tail(n_tail)) - tail_X;
                pearson += tail_diff * tail_diff / tail_X;
                n_bins += 1;
            }

            sum_chi += double(n_bins) * pearson;
            sum_df += double(n_bins - 1);

            if (return_envelope) {
                for (std::size_t p = 0; p < p_j; ++p) {
                    bin_probs.push_back(X(p));
                }
                if (n_tail > 0) {
                    bin_probs.push_back(tail_X);
                }
                bin_counts.push_back(int(n_bins));
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
