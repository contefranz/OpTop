#include <RcppArmadillo.h>
// [Rcpp::depends(RcppArmadillo)]

#ifdef _OPENMP
#include <omp.h>
#endif

// Blocked running-minimum kernel behind optop_make_partition(). For each
// vocabulary block it forms Theta_m %*% Phi_m[, block] for every model with
// one BLAS product, keeps the element-wise minimum across models in a
// J x block buffer, compares it with the per-document threshold tau_j, and
// writes the logical rare_mask block in place. The full J x W running-min
// matrix of the previous R implementation is never materialized, and no
// pmin() copies are allocated.
//
// The running minimum is initialized at the corpus baseline pi_glob(w): the
// harmonized rare set of the paper takes the union over the model grid AND
// the no-topics baseline, so a word is rare in document j iff
// min(pi_glob(w), min_K p^K_jw) < tau_j.
//
// Threading contract: within a block, the min update and the mask write are
// parallelized over columns, each column owned by exactly one thread, so the
// result is bit-identical for any n_threads (it is an exact comparison, with
// no floating-point reductions at all).

//' @keywords internal
// [[Rcpp::export]]
void optop_partition_fill_core(Rcpp::LogicalMatrix rare_mask,
                               const Rcpp::List& theta_list,
                               const Rcpp::List& phi_list,
                               const Rcpp::NumericVector& pi_glob,
                               const Rcpp::NumericVector& tau,
                               int block,
                               int n_threads)
{
    const arma::uword J = rare_mask.nrow();
    const arma::uword W = rare_mask.ncol();
    const int n_models = theta_list.size();

    if (phi_list.size() != n_models) {
        Rcpp::stop("theta_list and phi_list must have one entry per model");
    }
    if (static_cast<arma::uword>(pi_glob.size()) != W) {
        Rcpp::stop("pi_glob must have one entry per feature");
    }
    if (static_cast<arma::uword>(tau.size()) != J) {
        Rcpp::stop("tau must have one entry per document");
    }
    if (block < 1) {
        Rcpp::stop("block must be positive");
    }
    if (n_threads < 1) {
        n_threads = 1;
    }

    // no-copy views over the R matrices
    std::vector<arma::mat> thetas;
    std::vector<arma::mat> phis;
    thetas.reserve(n_models);
    phis.reserve(n_models);
    for (int m = 0; m < n_models; ++m) {
        Rcpp::NumericMatrix th = theta_list[m];
        Rcpp::NumericMatrix ph = phi_list[m];
        if (static_cast<arma::uword>(th.nrow()) != J) {
            Rcpp::stop("theta must have one row per document");
        }
        if (static_cast<arma::uword>(ph.ncol()) != W) {
            Rcpp::stop("phi must have one column per feature");
        }
        if (th.ncol() != ph.nrow()) {
            Rcpp::stop("theta and phi disagree on the number of topics");
        }
        thetas.emplace_back(th.begin(), th.nrow(), th.ncol(), false, true);
        phis.emplace_back(ph.begin(), ph.nrow(), ph.ncol(), false, true);
    }

    const double* tau_ptr = &tau[0];
    const double* pi_ptr = &pi_glob[0];
    int* mask_ptr = LOGICAL(rare_mask);
    arma::mat mn;

    for (arma::uword start = 0; start < W; start += static_cast<arma::uword>(block)) {
        Rcpp::checkUserInterrupt();
        const arma::uword end = std::min(start + static_cast<arma::uword>(block), W) - 1;
        const int wb = static_cast<int>(end - start + 1);

        // the baseline seeds the running minimum: column c holds the
        // constant pi_glob(start + c) before any model is folded in
        mn.set_size(J, wb);
        {
            double* dst = mn.memptr();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
            for (int c = 0; c < wb; ++c) {
                const double pw = pi_ptr[start + static_cast<std::size_t>(c)];
                const std::size_t off = static_cast<std::size_t>(c) * J;
                for (arma::uword j = 0; j < J; ++j) {
                    dst[off + j] = pw;
                }
            }
        }
        for (int m = 0; m < n_models; ++m) {
            const arma::mat I_m = thetas[m] * phis[m].cols(start, end);
            const double* src = I_m.memptr();
            double* dst = mn.memptr();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
            for (int c = 0; c < wb; ++c) {
                const std::size_t off = static_cast<std::size_t>(c) * J;
                for (arma::uword j = 0; j < J; ++j) {
                    if (src[off + j] < dst[off + j]) {
                        dst[off + j] = src[off + j];
                    }
                }
            }
        }

        const double* mn_ptr = mn.memptr();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
        for (int c = 0; c < wb; ++c) {
            const std::size_t src_off = static_cast<std::size_t>(c) * J;
            const std::size_t dst_off = (start + static_cast<std::size_t>(c)) * J;
            for (arma::uword j = 0; j < J; ++j) {
                mask_ptr[dst_off + j] = mn_ptr[src_off + j] < tau_ptr[j] ? 1 : 0;
            }
        }
    }
}

// Rare-mass kernel behind the Pearson min-bin inclusion rule. Given the
// harmonized mask, it accumulates for every model the fitted probability
// mass of the rare set, sum_{w in C*_j} p^m_jw, one vocabulary block at a
// time with the same blocked gemm as the fill kernel, and returns a
// J x n_models matrix. The caller scales by L_j and takes the minimum over
// models to evaluate min_K E^K_{j,min} >= c.
//
// Threading contract: within a block, documents are parallelized with each
// row owned by exactly one thread, and the per-document accumulation runs
// in ascending column order, so the result is bit-identical for any
// n_threads.

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix optop_partition_minmass_core(const Rcpp::LogicalMatrix& rare_mask,
                                                 const Rcpp::List& theta_list,
                                                 const Rcpp::List& phi_list,
                                                 int block,
                                                 int n_threads)
{
    const arma::uword J = rare_mask.nrow();
    const arma::uword W = rare_mask.ncol();
    const int n_models = theta_list.size();

    if (phi_list.size() != n_models) {
        Rcpp::stop("theta_list and phi_list must have one entry per model");
    }
    if (block < 1) {
        Rcpp::stop("block must be positive");
    }
    if (n_threads < 1) {
        n_threads = 1;
    }

    std::vector<arma::mat> thetas;
    std::vector<arma::mat> phis;
    thetas.reserve(n_models);
    phis.reserve(n_models);
    for (int m = 0; m < n_models; ++m) {
        Rcpp::NumericMatrix th = theta_list[m];
        Rcpp::NumericMatrix ph = phi_list[m];
        if (static_cast<arma::uword>(th.nrow()) != J) {
            Rcpp::stop("theta must have one row per document");
        }
        if (static_cast<arma::uword>(ph.ncol()) != W) {
            Rcpp::stop("phi must have one column per feature");
        }
        if (th.ncol() != ph.nrow()) {
            Rcpp::stop("theta and phi disagree on the number of topics");
        }
        thetas.emplace_back(th.begin(), th.nrow(), th.ncol(), false, true);
        phis.emplace_back(ph.begin(), ph.nrow(), ph.ncol(), false, true);
    }

    Rcpp::NumericMatrix out(static_cast<int>(J), n_models);
    double* o = REAL(out);
    const int* mask_ptr = LOGICAL(rare_mask);

    for (arma::uword start = 0; start < W; start += static_cast<arma::uword>(block)) {
        Rcpp::checkUserInterrupt();
        const arma::uword end = std::min(start + static_cast<arma::uword>(block), W) - 1;
        const int wb = static_cast<int>(end - start + 1);

        for (int m = 0; m < n_models; ++m) {
            const arma::mat I_m = thetas[m] * phis[m].cols(start, end);
            const double* ip = I_m.memptr();
            double* acc = o + static_cast<std::size_t>(m) * J;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
            for (arma::uword j = 0; j < J; ++j) {
                double s = 0.0;
                for (int c = 0; c < wb; ++c) {
                    if (mask_ptr[(start + static_cast<std::size_t>(c)) * J + j] == 1) {
                        s += ip[static_cast<std::size_t>(c) * J + j];
                    }
                }
                acc[j] += s;
            }
        }
    }
    return out;
}
