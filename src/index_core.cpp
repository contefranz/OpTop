#include <RcppArmadillo.h>
// [Rcpp::depends(RcppArmadillo)]

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Fused engine behind the OpTop discrepancy indices, format-2 generation.
// The document-level kernel consumes the sparse harmonized partition (the
// per-document non-rare word lists of partition_core.cpp) and merge-joins
// each list with the document's CSC column of the transposed counts, so the
// per-document cost is O(|NR_j| * K + nnz_j) with no J x W traversal, no
// packed mask, and no W x block gemm. The corpus and the partition cross
// the boundary as raw slot views (zero copy, no Armadillo sparse container,
// no shape limit).
//
// Per-family conventions (identical to the previous generation on the
// active support):
// * se: no eps flooring anywhere;
// * chisq: expected counts floored at eps element-wise on the ACTIVE
//   (non-rare) cells before any use;
// * deviance: contributions 2 * n * log(n / max(e, eps)) with the
//   0 * log(0) = 0 convention. At word level the FITTED deviance is the
//   Poisson unit deviance 2 * [n log(n/e) - (n - e)]: the linear correction
//   (unfloored e) makes every summand nonnegative, so the word-level index
//   never exceeds one; the correction is identically zero on the null side
//   and at document level (probability vectors on the support).
//
// Min-bin convention (format 2, one deliberate change): every collapsed-bin
// mass is the complement of a compensated non-rare sum, floored ONCE at eps
// where a denominator or a logarithm needs it: E_min = L_j * (1 - sum_NR p),
// B_min = L_j * (1 - sum_NR pi), N_min = the observed mass drained by the
// merge-join (exact). The previous generation floored the fitted chisq
// min-bin element-wise before summing (sum of max(e, eps) over rare cells),
// which cannot be reproduced without enumerating the rare cells; the
// difference is bounded by W * eps in absolute terms on a bin that enters
// the Pearson family only when its mass clears c, and the reference
// implementations follow the same single-floor convention.
//
// Threading contract (as everywhere in the package): each document (or
// word) is owned by exactly one thread and writes its own output slot;
// there are no cross-thread floating-point reductions, so results are
// bit-identical for any n_threads.
//
// Output layout of both kernels: a 6-column matrix
// [se, se_null, chisq, chisq_null, dev, dev_null], with columns of metrics
// (or sides) that were not requested left at zero.

namespace {

// Neumaier compensated accumulator (see partition_core.cpp)
struct KahanSum {
    double s = 0.0;
    double carry = 0.0;
    inline void add(const double v) {
        const double t = s + v;
        if (std::abs(s) >= std::abs(v)) {
            carry += (s - t) + v;
        } else {
            carry += (v - t) + s;
        }
        s = t;
    }
    inline double value() const { return s + carry; }
};

inline R_xlen_t off_at(const double* off, const R_xlen_t j)
{
    return static_cast<R_xlen_t>(off[j]);
}

}  // namespace

// Word-level kernel: the six accumulators of every word in the requested
// range, with the document dimension blocked INTERNALLY so that the fitted
// expectations exist only as a block x w_len gemm buffer (the previous
// generation received a full J x w_len dense block from R). Counts arrive
// as the raw CSC slots of the J x W matrix: each requested word is one
// contiguous column with ascending document indices, walked by a per-word
// cursor across the document blocks. Documents are accumulated in
// ascending order per word regardless of the schedule.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix optop_index_word_core(const arma::mat& theta,
                                          const arma::mat& phi_cols,
                                          const Rcpp::IntegerVector& N_p,
                                          const Rcpp::IntegerVector& N_i,
                                          const Rcpp::NumericVector& N_x,
                                          int w_start,
                                          int w_len,
                                          const Rcpp::NumericVector& L,
                                          const Rcpp::NumericVector& pi_w,
                                          double eps,
                                          bool do_model,
                                          bool do_null,
                                          bool do_se,
                                          bool do_chisq,
                                          bool do_dev,
                                          int n_threads)
{
    const R_xlen_t J = L.size();
    if (do_model) {
        if (theta.n_rows != static_cast<arma::uword>(J)) {
            Rcpp::stop("theta must have one row per document");
        }
        if (phi_cols.n_rows != theta.n_cols ||
            phi_cols.n_cols != static_cast<arma::uword>(w_len)) {
            Rcpp::stop("phi_cols does not match theta or the word range");
        }
    }
    if (pi_w.size() != w_len) {
        Rcpp::stop("pi_w must cover the requested word range");
    }
    if (static_cast<R_xlen_t>(N_p.size()) <
        static_cast<R_xlen_t>(w_start + w_len) + 1) {
        Rcpp::stop("the counts do not cover the requested word range");
    }
    if (n_threads < 1) n_threads = 1;

    Rcpp::NumericMatrix out(w_len, 6);
    double* o = REAL(out);
    const double* Lp = L.begin();
    const double* pip = pi_w.begin();
    const int* pp = N_p.begin();
    const int* ip = N_i.begin();
    const double* xp = N_x.begin();

    // per-word accumulators and sparse cursors, shared across blocks; each
    // word is touched by one thread per block and the blocks are serial
    std::vector<double> se(w_len, 0.0), se_n(w_len, 0.0);
    std::vector<double> x2(w_len, 0.0), x2_n(w_len, 0.0);
    std::vector<double> dv(w_len, 0.0), dv_n(w_len, 0.0);
    std::vector<double> dv_lin(w_len, 0.0);
    std::vector<int> cursor(w_len);
    for (int c = 0; c < w_len; ++c) {
        cursor[c] = pp[w_start + c];
    }

    const R_xlen_t block_docs = 8192;
    arma::mat E_sub;

    for (R_xlen_t j0 = 0; j0 < J; j0 += block_docs) {
        Rcpp::checkUserInterrupt();
        const R_xlen_t j1 = std::min<R_xlen_t>(j0 + block_docs, J);

        if (do_model) {
            // block x w_len fitted probabilities in one gemm
            E_sub = theta.rows(static_cast<arma::uword>(j0),
                               static_cast<arma::uword>(j1 - 1)) * phi_cols;
        }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) num_threads(n_threads)
#endif
        for (int c = 0; c < w_len; ++c) {
            const int col_end = pp[w_start + c + 1];
            int cur = cursor[c];
            const double pw = pip[c];
            const double* e_col = do_model ?
                E_sub.colptr(static_cast<arma::uword>(c)) : nullptr;

            double a_se = 0.0, a_se_n = 0.0;
            double a_x2 = 0.0, a_x2_n = 0.0;
            double a_dv = 0.0, a_dv_n = 0.0;
            double a_dv_lin = 0.0;

            for (R_xlen_t j = j0; j < j1; ++j) {
                double n = 0.0;
                if (cur < col_end && static_cast<R_xlen_t>(ip[cur]) == j) {
                    n = xp[cur];
                    ++cur;
                }
                if (do_model) {
                    const double e = Lp[j] * e_col[j - j0];
                    if (do_se) {
                        const double d = n - e;
                        a_se += d * d;
                    }
                    if (do_chisq) {
                        const double ef = e < eps ? eps : e;
                        const double d = n - ef;
                        a_x2 += d * d / ef;
                    }
                    if (do_dev) {
                        // Poisson unit deviance: the linear term uses the
                        // unfloored expectation and runs over every document
                        a_dv_lin += n - e;
                        if (n > 0.0) {
                            const double ef = e < eps ? eps : e;
                            a_dv += n * (std::log(n) - std::log(ef));
                        }
                    }
                }
                if (do_null) {
                    const double b = Lp[j] * pw;
                    if (do_se) {
                        const double d = n - b;
                        a_se_n += d * d;
                    }
                    if (do_chisq) {
                        const double bf = b < eps ? eps : b;
                        const double d = n - bf;
                        a_x2_n += d * d / bf;
                    }
                    if (do_dev && n > 0.0) {
                        const double bf = b < eps ? eps : b;
                        a_dv_n += n * (std::log(n) - std::log(bf));
                    }
                }
            }

            se[c] += a_se;
            se_n[c] += a_se_n;
            x2[c] += a_x2;
            x2_n[c] += a_x2_n;
            dv[c] += a_dv;
            dv_n[c] += a_dv_n;
            dv_lin[c] += a_dv_lin;
            cursor[c] = cur;
        }
    }

    for (int c = 0; c < w_len; ++c) {
        o[c] = se[c];
        o[c + static_cast<std::size_t>(w_len)] = se_n[c];
        o[c + 2 * static_cast<std::size_t>(w_len)] = x2[c];
        o[c + 3 * static_cast<std::size_t>(w_len)] = x2_n[c];
        o[c + 4 * static_cast<std::size_t>(w_len)] = 2.0 * (dv[c] - dv_lin[c]);
        o[c + 5 * static_cast<std::size_t>(w_len)] = 2.0 * dv_n[c];
    }
    return out;
}

// Document-level kernel on the harmonized support {non-rare terms} U {min}.
// One merge-join per document: the sorted non-rare list against the
// document's CSC column of the transposed counts (word indices ascending).
// Non-rare cells accumulate their per-family terms DIRECTLY, in ascending
// word order (a zero-count non-rare cell contributes its zero-count value);
// observed cells outside the list drain into the exact observed min-bin
// mass N_min; the fitted and baseline min-bin masses follow by compensated
// complement. Per-document cost O(|NR_j| * K + nnz_j).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix optop_index_doc_core(const arma::mat& theta,
                                         const arma::mat& phi,
                                         const Rcpp::IntegerVector& Nt_p,
                                         const Rcpp::IntegerVector& Nt_i,
                                         const Rcpp::NumericVector& Nt_x,
                                         const Rcpp::NumericVector& nr_off,
                                         const Rcpp::IntegerVector& nr_words,
                                         const Rcpp::NumericVector& L,
                                         const Rcpp::NumericVector& pi_row,
                                         const Rcpp::LogicalVector& chisq_min_ok,
                                         double eps,
                                         bool do_model,
                                         bool do_null,
                                         bool do_se,
                                         bool do_chisq,
                                         bool do_dev,
                                         int n_threads)
{
    const R_xlen_t J = L.size();
    const arma::uword K = do_model ? theta.n_cols : 0;
    if (Nt_p.size() != J + 1) {
        Rcpp::stop("the transposed counts do not match the document set");
    }
    if (nr_off.size() != J + 1) {
        Rcpp::stop("the partition does not match the document set");
    }
    if (chisq_min_ok.size() != J) {
        Rcpp::stop("chisq_min_ok must have one entry per document");
    }
    if (do_model) {
        if (theta.n_rows != static_cast<arma::uword>(J)) {
            Rcpp::stop("theta must have one row per document");
        }
        if (phi.n_rows != K) {
            Rcpp::stop("theta and phi disagree on the number of topics");
        }
    }
    if (n_threads < 1) n_threads = 1;

    Rcpp::NumericMatrix out(J, 6);
    double* o = REAL(out);
    const double* Lp = L.begin();
    const double* pip = pi_row.begin();
    const int* min_ok = LOGICAL(chisq_min_ok);
    const int* pp = Nt_p.begin();
    const int* ipv = Nt_i.begin();
    const double* xpv = Nt_x.begin();
    const double* offp = nr_off.begin();
    const int* wp = nr_words.begin();

    // serial macro-chunks keep the loop interruptible; documents inside a
    // chunk run in parallel, each writing its own row
    const R_xlen_t chunk = 65536;
    for (R_xlen_t c0 = 0; c0 < J; c0 += chunk) {
        Rcpp::checkUserInterrupt();
        const R_xlen_t c1 = std::min<R_xlen_t>(c0 + chunk, J);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
        {
            std::vector<double> th(K > 0 ? K : 1);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
            for (R_xlen_t j = c0; j < c1; ++j) {
                const double Lj = Lp[j];
                if (do_model) {
                    for (arma::uword k = 0; k < K; ++k) {
                        th[k] = theta.at(static_cast<arma::uword>(j), k);
                    }
                }

                double se_nr = 0.0, seN_nr = 0.0;
                double x2_nr = 0.0, x2N_nr = 0.0;
                double dv_nr = 0.0, dvN_nr = 0.0;
                double N_min = 0.0;
                KahanSum psum;
                KahanSum pisum;

                R_xlen_t a = off_at(offp, j);
                const R_xlen_t a1 = off_at(offp, j + 1);
                int b = pp[j];
                const int b1 = pp[j + 1];

                while (a < a1) {
                    const int w = wp[a];
                    // rare observed cells before the next non-rare word
                    while (b < b1 && ipv[b] < w) {
                        N_min += xpv[b];
                        ++b;
                    }
                    double n = 0.0;
                    if (b < b1 && ipv[b] == w) {
                        n = xpv[b];
                        ++b;
                    }

                    if (do_model) {
                        double p = 0.0;
                        const double* ph = phi.colptr(
                            static_cast<arma::uword>(w));
                        for (arma::uword k = 0; k < K; ++k) {
                            p += th[k] * ph[k];
                        }
                        psum.add(p);
                        const double e = Lj * p;
                        if (do_se) {
                            const double d = n - e;
                            se_nr += d * d;
                        }
                        if (do_chisq) {
                            const double ef = e < eps ? eps : e;
                            const double d = n - ef;
                            x2_nr += d * d / ef;
                        }
                        if (do_dev && n > 0.0) {
                            const double ef = e < eps ? eps : e;
                            dv_nr += 2.0 * n * std::log(n / ef);
                        }
                    }
                    if (do_null) {
                        const double piw = pip[w];
                        pisum.add(piw);
                        const double bb = Lj * piw;
                        if (do_se) {
                            const double d = n - bb;
                            seN_nr += d * d;
                        }
                        if (do_chisq) {
                            const double bf = bb < eps ? eps : bb;
                            const double d = n - bf;
                            x2N_nr += d * d / bf;
                        }
                        if (do_dev && n > 0.0) {
                            const double bf = bb < eps ? eps : bb;
                            dvN_nr += 2.0 * n * std::log(n / bf);
                        }
                    }
                    ++a;
                }
                // remaining observed cells are rare
                while (b < b1) {
                    N_min += xpv[b];
                    ++b;
                }

                // min-bin masses by compensated complement: exact when the
                // non-rare list is empty (structural collapse), absolute
                // error O(L_j * eps_machine) otherwise
                const double E_min = do_model ?
                    Lj * std::max(0.0, 1.0 - psum.value()) : 0.0;
                const double B_min = do_null ?
                    Lj * std::max(0.0, 1.0 - pisum.value()) : 0.0;

                if (do_model) {
                    if (do_se) {
                        const double d = N_min - E_min;
                        o[j] = se_nr + d * d;
                    }
                    if (do_chisq) {
                        // Pearson inclusion rule: the min bin enters only
                        // when the grid-wide flag holds; excluded bins drop
                        // from BOTH sides
                        if (min_ok[j] == 1) {
                            const double em = E_min < eps ? eps : E_min;
                            const double d = N_min - em;
                            o[j + 2 * static_cast<std::size_t>(J)] =
                                x2_nr + d * d / em;
                        } else {
                            o[j + 2 * static_cast<std::size_t>(J)] = x2_nr;
                        }
                    }
                    if (do_dev) {
                        const double em = E_min < eps ? eps : E_min;
                        const double nm = N_min < eps ? eps : N_min;
                        const double dev_min = N_min == 0.0 ?
                            0.0 : 2.0 * N_min * std::log(nm / em);
                        o[j + 4 * static_cast<std::size_t>(J)] =
                            dv_nr + dev_min;
                    }
                }
                if (do_null) {
                    if (do_se) {
                        const double d = N_min - B_min;
                        o[j + static_cast<std::size_t>(J)] = seN_nr + d * d;
                    }
                    if (do_chisq) {
                        if (min_ok[j] == 1) {
                            const double bm = B_min < eps ? eps : B_min;
                            const double d = N_min - bm;
                            o[j + 3 * static_cast<std::size_t>(J)] =
                                x2N_nr + d * d / bm;
                        } else {
                            o[j + 3 * static_cast<std::size_t>(J)] = x2N_nr;
                        }
                    }
                    if (do_dev) {
                        const double bm = B_min < eps ? eps : B_min;
                        const double nm = N_min < eps ? eps : N_min;
                        const double dev_min = N_min == 0.0 ?
                            0.0 : 2.0 * N_min * std::log(nm / bm);
                        o[j + 5 * static_cast<std::size_t>(J)] =
                            dvN_nr + dev_min;
                    }
                }
            }
        }
    }

    return out;
}
