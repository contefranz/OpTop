#include <RcppArmadillo.h>
// [Rcpp::depends(RcppArmadillo)]

#include <cmath>
#include <cstdint>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Fused engine behind the OpTop discrepancy indices. The R implementation
// evaluated each metric with a chain of J x W temporaries per model
// (difference matrices, logical-to-double mask coercions, outer() baselines,
// masked rowSums); at the scales of the efficiency study those temporaries,
// not the BLAS product, dominate the runtime and the memory profile. The
// kernels below reproduce the exact conventions of the R code (see
// R/discrepancies.R and the naive references in tests/testthat) in a single
// traversal per document or word, computing every requested metric and, when
// asked, the model-independent baseline side at the same time:
// * se: no eps flooring anywhere;
// * chisq: expected counts floored at eps element-wise before any use; the
//   collapsed min bin uses the sum of the floored elements, floored again at
//   eps; the baseline min bin is L_j * sum(pi[rare]) from the unfloored pi,
//   then floored;
// * deviance: contributions 2 * n * log(n / max(e, eps)) with the
//   0 * log(0) = 0 convention; the min bin uses the sum of the unfloored
//   expected counts, floored at eps. At word level the FITTED deviance is
//   the Poisson unit deviance 2 * [n log(n/e) - (n - e)]: the linear
//   correction (unfloored e) makes every summand nonnegative, so the
//   word-level index never exceeds one. The correction is identically zero
//   on the null side because sum_j B_jw = sum_j N_jw for the in-sample
//   corpus baseline, and at document level because the binned fitted and
//   empirical vectors are both probability vectors on the support.
//
// Threading contract (as everywhere in the package): each document (or word)
// is owned by exactly one thread and writes its own output slot; there are
// no cross-thread floating-point reductions, so results are bit-identical
// for any n_threads.
//
// Output layout of both kernels: a 6-column matrix
// [se, se_null, chisq, chisq_null, dev, dev_null], with columns of metrics
// (or sides) that were not requested left at zero.

namespace {

inline bool bit_get(const unsigned char* bits, const std::size_t idx)
{
    return (bits[idx >> 3] >> (idx & 7)) & 1u;
}

}  // namespace

// Pack a J x W logical mask into bits indexed by (j * W + w), so that a
// document's mask is a contiguous bit range. One pass; ~J*W/8 bytes.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::RawVector optop_pack_mask_core(const Rcpp::LogicalMatrix& mask)
{
    const std::size_t J = mask.nrow();
    const std::size_t W = mask.ncol();
    Rcpp::RawVector out((J * W + 7) / 8);
    unsigned char* bits = RAW(out);
    std::memset(bits, 0, out.size());

    const int* m = LOGICAL(mask);
    for (std::size_t w = 0; w < W; ++w) {
        const std::size_t col = w * J;
        for (std::size_t j = 0; j < J; ++j) {
            if (m[col + j] == 1) {
                const std::size_t idx = j * W + w;
                bits[idx >> 3] |= static_cast<unsigned char>(1u << (idx & 7));
            }
        }
    }
    return out;
}

// Word-level kernel: one word (one dense column of E, one sparse column of
// N) per iteration. E_block is the J x wb fitted-count block computed in R
// by BLAS, unfloored; the per-metric floors are applied here. The harmonized
// partition plays no role at word level.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix optop_index_word_core(const arma::mat& E_block,
                                          const arma::sp_mat& N,
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
    const std::size_t J = N.n_rows;
    if (do_model && (E_block.n_rows != J ||
                     E_block.n_cols != static_cast<arma::uword>(w_len))) {
        Rcpp::stop("E_block does not match the requested word range");
    }
    if (static_cast<std::size_t>(L.size()) != J) {
        Rcpp::stop("L must have one entry per document");
    }
    if (pi_w.size() != w_len) {
        Rcpp::stop("pi_w must cover the requested word range");
    }
    if (n_threads < 1) n_threads = 1;

    Rcpp::NumericMatrix out(w_len, 6);
    double* o = REAL(out);
    const double* Lp = &L[0];
    const double* pip = &pi_w[0];

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
    {
        // thread-local dense scatter of one sparse column
        std::vector<double> n_dense(J, 0.0);
        std::vector<arma::uword> touched;
        touched.reserve(256);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
        for (int c = 0; c < w_len; ++c) {
            const arma::uword w = static_cast<arma::uword>(w_start + c);
            for (arma::sp_mat::const_col_iterator it = N.begin_col(w);
                 it != N.end_col(w); ++it) {
                n_dense[it.row()] = *it;
                touched.push_back(it.row());
            }

            const double pw = pip[c];
            const double* e_col = do_model ? E_block.colptr(c) : nullptr;

            double se = 0.0, se_n = 0.0;
            double x2 = 0.0, x2_n = 0.0;
            double dv = 0.0, dv_n = 0.0;
            double dv_lin = 0.0;

            for (std::size_t j = 0; j < J; ++j) {
                const double n = n_dense[j];
                if (do_model) {
                    const double e = e_col[j];
                    if (do_se) {
                        const double d = n - e;
                        se += d * d;
                    }
                    if (do_chisq) {
                        const double ef = e < eps ? eps : e;
                        const double d = n - ef;
                        x2 += d * d / ef;
                    }
                    if (do_dev) {
                        // Poisson unit deviance: the linear term uses the
                        // unfloored expectation and runs over every document
                        dv_lin += n - e;
                        if (n > 0.0) {
                            const double ef = e < eps ? eps : e;
                            dv += n * (std::log(n) - std::log(ef));
                        }
                    }
                }
                if (do_null) {
                    const double b = Lp[j] * pw;
                    if (do_se) {
                        const double d = n - b;
                        se_n += d * d;
                    }
                    if (do_chisq) {
                        const double bf = b < eps ? eps : b;
                        const double d = n - bf;
                        x2_n += d * d / bf;
                    }
                    if (do_dev && n > 0.0) {
                        const double bf = b < eps ? eps : b;
                        dv_n += n * (std::log(n) - std::log(bf));
                    }
                }
            }

            o[c] = se;
            o[c + static_cast<std::size_t>(w_len)] = se_n;
            o[c + 2 * static_cast<std::size_t>(w_len)] = x2;
            o[c + 3 * static_cast<std::size_t>(w_len)] = x2_n;
            o[c + 4 * static_cast<std::size_t>(w_len)] = 2.0 * (dv - dv_lin);
            o[c + 5 * static_cast<std::size_t>(w_len)] = 2.0 * dv_n;

            for (const arma::uword j : touched) {
                n_dense[j] = 0.0;
            }
            touched.clear();
        }
    }
    return out;
}

// Document-level kernel: one document per iteration, on the harmonized
// support {non-rare terms} U {min}. The fitted probabilities of the block
// are produced by one gemm inside the kernel (tww is phi transposed, W x K),
// so a document is a contiguous column; counts arrive as the transposed
// sparse dfm (documents are CSC columns) and the rare mask as packed bits.
// The dense pass accumulates the zero-count value of every term and the
// collapsed-bin sums; the sparse pass replaces the zero-count value with the
// observed-count value at the nonzeros.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix optop_index_doc_core(const arma::mat& tww,
                                         const arma::mat& theta_blk,
                                         const arma::sp_mat& N_t,
                                         int doc_start,
                                         const Rcpp::RawVector& mask_bits,
                                         const Rcpp::NumericVector& L_blk,
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
    const std::size_t W = N_t.n_rows;
    const int b = L_blk.size();
    if (chisq_min_ok.size() != b) {
        Rcpp::stop("chisq_min_ok must have one entry per document of the block");
    }
    if (do_model) {
        if (tww.n_rows != W) {
            Rcpp::stop("tww must have one row per feature");
        }
        if (theta_blk.n_rows != static_cast<arma::uword>(b) ||
            theta_blk.n_cols != tww.n_cols) {
            Rcpp::stop("theta_blk does not match tww or the block length");
        }
    }
    if (static_cast<std::size_t>(pi_row.size()) != W) {
        Rcpp::stop("pi_row must have one entry per feature");
    }
    if (n_threads < 1) n_threads = 1;

    // one gemm per block: fitted probabilities, one column per document
    arma::mat I_t;
    if (do_model) {
        I_t = tww * theta_blk.t();  // W x b
    }

    Rcpp::NumericMatrix out(b, 6);
    double* o = REAL(out);
    const double* Lp = &L_blk[0];
    const double* pip = &pi_row[0];
    const int* min_ok = LOGICAL(chisq_min_ok);
    const unsigned char* bits = RAW(mask_bits);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
    for (int jj = 0; jj < b; ++jj) {
        const std::size_t j_global = static_cast<std::size_t>(doc_start + jj);
        const std::size_t bit0 = j_global * W;
        const double Lj = Lp[jj];
        const double* i_col = do_model ? I_t.colptr(jj) : nullptr;

        // dense pass: zero-count values over the full support, plus the
        // collapsed-bin sums (ascending term order, as in the R rowSums)
        double se_full = 0.0, se_rare = 0.0, E_min = 0.0;
        double x2_full = 0.0, x2_rare = 0.0, E_min_f = 0.0;
        double seN_full = 0.0, seN_rare = 0.0;
        double x2N_full = 0.0, x2N_rare = 0.0;
        double pi_rare = 0.0;

        for (std::size_t w = 0; w < W; ++w) {
            const bool rare = bit_get(bits, bit0 + w);
            if (do_model) {
                const double e = i_col[w] * Lj;
                if (do_se || do_dev) {
                    if (do_se) {
                        const double v = e * e;
                        se_full += v;
                        if (rare) se_rare += v;
                    }
                    if (rare) E_min += e;
                }
                if (do_chisq) {
                    const double ef = e < eps ? eps : e;
                    x2_full += ef;
                    if (rare) {
                        x2_rare += ef;
                        E_min_f += ef;
                    }
                }
            }
            if (do_null) {
                if (rare) pi_rare += pip[w];
                const double bb = Lj * pip[w];
                if (do_se) {
                    const double v = bb * bb;
                    seN_full += v;
                    if (rare) seN_rare += v;
                }
                if (do_chisq) {
                    const double bf = bb < eps ? eps : bb;
                    x2N_full += bf;
                    if (rare) x2N_rare += bf;
                }
            }
        }

        // sparse pass: replace the zero-count value with the observed one at
        // the nonzeros, and collect the observed min-bin mass
        double N_min = 0.0;
        double dv_full = 0.0, dv_rare = 0.0;
        double dvN_full = 0.0, dvN_rare = 0.0;

        for (arma::sp_mat::const_col_iterator it = N_t.begin_col(j_global);
             it != N_t.end_col(j_global); ++it) {
            const std::size_t w = it.row();
            const double n = *it;
            const bool rare = bit_get(bits, bit0 + w);
            if (rare) N_min += n;

            if (do_model) {
                const double e = i_col[w] * Lj;
                if (do_se) {
                    const double d = n - e;
                    const double corr = d * d - e * e;
                    se_full += corr;
                    if (rare) se_rare += corr;
                }
                if (do_chisq) {
                    const double ef = e < eps ? eps : e;
                    const double d = n - ef;
                    const double corr = d * d / ef - ef;
                    x2_full += corr;
                    if (rare) x2_rare += corr;
                }
                if (do_dev && n > 0.0) {
                    const double ef = e < eps ? eps : e;
                    const double contrib = 2.0 * n * std::log(n / ef);
                    dv_full += contrib;
                    if (rare) dv_rare += contrib;
                }
            }
            if (do_null) {
                const double bb = Lj * pip[w];
                if (do_se) {
                    const double d = n - bb;
                    const double corr = d * d - bb * bb;
                    seN_full += corr;
                    if (rare) seN_rare += corr;
                }
                if (do_chisq) {
                    const double bf = bb < eps ? eps : bb;
                    const double d = n - bf;
                    const double corr = d * d / bf - bf;
                    x2N_full += corr;
                    if (rare) x2N_rare += corr;
                }
                if (do_dev && n > 0.0) {
                    const double bf = bb < eps ? eps : bb;
                    const double contrib = 2.0 * n * std::log(n / bf);
                    dvN_full += contrib;
                    if (rare) dvN_rare += contrib;
                }
            }
        }

        const double B_min_raw = Lj * pi_rare;

        if (do_model) {
            if (do_se) {
                const double d = N_min - E_min;
                o[jj] = se_full - se_rare + d * d;
            }
            if (do_chisq) {
                // Pearson inclusion rule: the min bin enters only when the
                // grid-wide flag holds; excluded bins drop from BOTH sides
                if (min_ok[jj] == 1) {
                    const double em = E_min_f < eps ? eps : E_min_f;
                    const double d = N_min - em;
                    o[jj + 2 * static_cast<std::size_t>(b)] =
                        x2_full - x2_rare + d * d / em;
                } else {
                    o[jj + 2 * static_cast<std::size_t>(b)] =
                        x2_full - x2_rare;
                }
            }
            if (do_dev) {
                const double em = E_min < eps ? eps : E_min;
                const double nm = N_min < eps ? eps : N_min;
                const double dev_min =
                    N_min == 0.0 ? 0.0 : 2.0 * N_min * std::log(nm / em);
                o[jj + 4 * static_cast<std::size_t>(b)] =
                    dv_full - dv_rare + dev_min;
            }
        }
        if (do_null) {
            if (do_se) {
                const double d = N_min - B_min_raw;
                o[jj + static_cast<std::size_t>(b)] =
                    seN_full - seN_rare + d * d;
            }
            if (do_chisq) {
                if (min_ok[jj] == 1) {
                    const double bm = B_min_raw < eps ? eps : B_min_raw;
                    const double d = N_min - bm;
                    o[jj + 3 * static_cast<std::size_t>(b)] =
                        x2N_full - x2N_rare + d * d / bm;
                } else {
                    o[jj + 3 * static_cast<std::size_t>(b)] =
                        x2N_full - x2N_rare;
                }
            }
            if (do_dev) {
                const double bm = B_min_raw < eps ? eps : B_min_raw;
                const double nm = N_min < eps ? eps : N_min;
                const double dev_min =
                    N_min == 0.0 ? 0.0 : 2.0 * N_min * std::log(nm / bm);
                o[jj + 5 * static_cast<std::size_t>(b)] =
                    dvN_full - dvN_rare + dev_min;
            }
        }
    }

    return out;
}
