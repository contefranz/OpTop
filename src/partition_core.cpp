#include <RcppArmadillo.h>
// [Rcpp::depends(RcppArmadillo)]

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Sparse harmonized partition (format 2). The partition stores the
// COMPLEMENT of the rare set: the per-document non-rare word lists, whose
// total size is bounded by sum_j L_j / c because a document can hold at most
// L_j / c words with p_jw >= tau_j = c / L_j (probabilities sum to one). The
// dense J x W rare mask of earlier versions is never materialized.
//
// The construction is exact and rests on the augmented-union rule of the
// paper (K0 = K U {null}): a cell is non-rare only if EVERY member of the
// augmented grid clears tau_j, in particular the null baseline, so
// pi_glob(w) >= tau_j is a necessary condition. Sorting the vocabulary once
// by pi_glob descending, each document's candidate set is the prefix with
// pi_glob >= tau_j (at most L_j / c words), and only candidates need the
// per-model check p^K_jw >= tau_j. Total cost O((tokens / c) * sum_K K)
// instead of the former O(J * W * sum_K K) blocked products.
//
// The per-model pass is a separate kernel so that models can be
// materialized one at a time (lazy loading): pass 1 ANDs the per-model
// keep flags over the shared candidate prefixes; after compaction, pass 2
// re-evaluates each model on the final non-rare lists to accumulate
// sum_{w in NR_j} p^K_jw, from which the min-bin masses follow by
// complement: E^K_{j,min} = L_j * (1 - sum_{NR_j} p^K_jw). Complement sums
// use Neumaier compensation, so the absolute error is O(eps_machine) on the
// probability scale; a document with an empty non-rare list has sums of
// exactly zero and therefore exact min-bin masses (the structural-collapse
// case of the null-discrepancy floor stays exact).
//
// Threading contract (package-wide): each document is owned by exactly one
// thread and writes its own slots; no cross-thread floating-point
// reductions, so results are bit-identical for any n_threads.

namespace {

// Neumaier compensated accumulator: absolute error O(eps_machine) of the
// total instead of O(n * eps_machine)
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

// Candidate prefix lengths: pi_sorted holds pi_glob sorted DESCENDING and
// prefix_len[j] = #{t : pi_sorted[t] >= tau_j}, located by binary search.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::IntegerVector optop_partition_candidates_core(
        const Rcpp::NumericVector& pi_sorted,
        const Rcpp::NumericVector& tau,
        int n_threads)
{
    const R_xlen_t J = tau.size();
    const R_xlen_t W = pi_sorted.size();
    Rcpp::IntegerVector prefix_len(J);
    const double* ps = pi_sorted.begin();
    const double* tp = tau.begin();
    int* out = prefix_len.begin();
    if (n_threads < 1) n_threads = 1;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (R_xlen_t j = 0; j < J; ++j) {
        // first position with pi_sorted < tau_j on the descending array
        const double* it = std::lower_bound(
            ps, ps + W, tp[j],
            [](const double a, const double b) { return a >= b; });
        out[j] = static_cast<int>(it - ps);
    }
    return prefix_len;
}

// One model's pass over the shared candidate prefixes: keep[t] survives only
// while every model seen so far satisfies p^K_jw >= tau_j at the candidate
// cell. keep is a raw vector updated in place (the internal caller owns it;
// candidate t of document j is the word order[t - cand_off[j]]).
//' @keywords internal
// [[Rcpp::export]]
void optop_partition_pass_core(const arma::mat& theta,
                               const arma::mat& phi,
                               const Rcpp::IntegerVector& order,
                               const Rcpp::NumericVector& cand_off,
                               Rcpp::RawVector keep,
                               const Rcpp::NumericVector& tau,
                               int n_threads)
{
    const R_xlen_t J = tau.size();
    const arma::uword K = theta.n_cols;
    if (phi.n_rows != K) {
        Rcpp::stop("theta and phi disagree on the number of topics");
    }
    if (theta.n_rows != static_cast<arma::uword>(J)) {
        Rcpp::stop("theta must have one row per document");
    }
    if (cand_off.size() != J + 1) {
        Rcpp::stop("cand_off must have J + 1 entries");
    }
    const int* ord = order.begin();
    const double* offp = cand_off.begin();
    const double* tp = tau.begin();
    unsigned char* kp = RAW(keep);
    if (n_threads < 1) n_threads = 1;

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
    {
        std::vector<double> th(K);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
        for (R_xlen_t j = 0; j < J; ++j) {
            const R_xlen_t o0 = off_at(offp, j);
            const R_xlen_t o1 = off_at(offp, j + 1);
            if (o0 == o1) continue;
            for (arma::uword k = 0; k < K; ++k) {
                th[k] = theta.at(static_cast<arma::uword>(j), k);
            }
            const double tau_j = tp[j];
            for (R_xlen_t t = o0; t < o1; ++t) {
                if (!kp[t]) continue;
                const double* ph = phi.colptr(
                    static_cast<arma::uword>(ord[t - o0]));
                double p = 0.0;
                for (arma::uword k = 0; k < K; ++k) {
                    p += th[k] * ph[k];
                }
                if (p < tau_j) kp[t] = 0;
            }
        }
    }
}

// Compact the surviving candidates into the ragged non-rare structure:
// offsets (double, exact to 2^53) and word indices, sorted ascending per
// document so the index kernels can merge-join against CSC columns.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List optop_partition_compact_core(const Rcpp::NumericVector& cand_off,
                                        const Rcpp::RawVector& keep,
                                        const Rcpp::IntegerVector& order,
                                        int n_threads)
{
    const R_xlen_t J = cand_off.size() - 1;
    const double* offp = cand_off.begin();
    const unsigned char* kp = RAW(keep);
    const int* ord = order.begin();
    if (n_threads < 1) n_threads = 1;

    // pass 1: per-document survivor counts
    std::vector<R_xlen_t> counts(J, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (R_xlen_t j = 0; j < J; ++j) {
        R_xlen_t n = 0;
        for (R_xlen_t t = off_at(offp, j); t < off_at(offp, j + 1); ++t) {
            n += kp[t] ? 1 : 0;
        }
        counts[j] = n;
    }

    Rcpp::NumericVector nr_off(J + 1);
    double* no = nr_off.begin();
    no[0] = 0.0;
    for (R_xlen_t j = 0; j < J; ++j) {
        no[j + 1] = no[j] + static_cast<double>(counts[j]);
    }

    const R_xlen_t total = static_cast<R_xlen_t>(no[J]);
    Rcpp::IntegerVector nr_words(total);
    int* wp = nr_words.begin();

    // pass 2: fill and sort ascending per document
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(n_threads)
#endif
    for (R_xlen_t j = 0; j < J; ++j) {
        R_xlen_t dst = off_at(no, j);
        const R_xlen_t first = dst;
        const R_xlen_t o0 = off_at(offp, j);
        const R_xlen_t o1 = off_at(offp, j + 1);
        for (R_xlen_t t = o0; t < o1; ++t) {
            if (kp[t]) {
                wp[dst++] = ord[t - o0];
            }
        }
        std::sort(wp + first, wp + dst);
    }

    return Rcpp::List::create(Rcpp::Named("offsets") = nr_off,
                              Rcpp::Named("words") = nr_words);
}

// Per-document compensated sum of one model's fitted probabilities over the
// final non-rare list: sum_{w in NR_j} p^K_jw. The caller derives the
// min-bin fitted mass by complement, E^K_{j,min} = L_j * (1 - sum).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector optop_partition_sums_core(const arma::mat& theta,
                                              const arma::mat& phi,
                                              const Rcpp::NumericVector& nr_off,
                                              const Rcpp::IntegerVector& nr_words,
                                              int n_threads)
{
    const R_xlen_t J = nr_off.size() - 1;
    const arma::uword K = theta.n_cols;
    if (phi.n_rows != K) {
        Rcpp::stop("theta and phi disagree on the number of topics");
    }
    if (theta.n_rows != static_cast<arma::uword>(J)) {
        Rcpp::stop("theta must have one row per document");
    }
    Rcpp::NumericVector sums(J);
    const double* offp = nr_off.begin();
    const int* wp = nr_words.begin();
    double* out = sums.begin();
    if (n_threads < 1) n_threads = 1;

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
    {
        std::vector<double> th(K);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
        for (R_xlen_t j = 0; j < J; ++j) {
            const R_xlen_t o0 = off_at(offp, j);
            const R_xlen_t o1 = off_at(offp, j + 1);
            if (o0 == o1) {
                out[j] = 0.0;
                continue;
            }
            for (arma::uword k = 0; k < K; ++k) {
                th[k] = theta.at(static_cast<arma::uword>(j), k);
            }
            KahanSum acc;
            for (R_xlen_t t = o0; t < o1; ++t) {
                const double* ph = phi.colptr(static_cast<arma::uword>(wp[t]));
                double p = 0.0;
                for (arma::uword k = 0; k < K; ++k) {
                    p += th[k] * ph[k];
                }
                acc.add(p);
            }
            out[j] = acc.value();
        }
    }
    return sums;
}

// Per-document compensated sum of the baseline probabilities over the
// non-rare list: sum_{w in NR_j} pi_glob(w), for B_{j,min} by complement.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector optop_partition_pisum_core(const Rcpp::NumericVector& nr_off,
                                               const Rcpp::IntegerVector& nr_words,
                                               const Rcpp::NumericVector& pi_glob,
                                               int n_threads)
{
    const R_xlen_t J = nr_off.size() - 1;
    Rcpp::NumericVector sums(J);
    const double* offp = nr_off.begin();
    const int* wp = nr_words.begin();
    const double* pip = pi_glob.begin();
    double* out = sums.begin();
    if (n_threads < 1) n_threads = 1;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (R_xlen_t j = 0; j < J; ++j) {
        KahanSum acc;
        for (R_xlen_t t = off_at(offp, j); t < off_at(offp, j + 1); ++t) {
            acc.add(pip[wp[t]]);
        }
        out[j] = acc.value();
    }
    return sums;
}

// Per-document observed mass on the non-rare list, by merge-join of the
// sorted word list with the document's CSC column of the transposed counts
// (word indices ascending): sum_{w in NR_j} N_jw. The observed min-bin mass
// follows as N_{j,min} = L_j - sum.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector optop_partition_obsmass_core(const Rcpp::IntegerVector& Nt_p,
                                                 const Rcpp::IntegerVector& Nt_i,
                                                 const Rcpp::NumericVector& Nt_x,
                                                 const Rcpp::NumericVector& nr_off,
                                                 const Rcpp::IntegerVector& nr_words,
                                                 int n_threads)
{
    const R_xlen_t J = nr_off.size() - 1;
    if (Nt_p.size() != J + 1) {
        Rcpp::stop("the transposed counts do not match the partition");
    }
    Rcpp::NumericVector sums(J);
    const int* pp = Nt_p.begin();
    const int* ip = Nt_i.begin();
    const double* xp = Nt_x.begin();
    const double* offp = nr_off.begin();
    const int* wp = nr_words.begin();
    double* out = sums.begin();
    if (n_threads < 1) n_threads = 1;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (R_xlen_t j = 0; j < J; ++j) {
        R_xlen_t a = off_at(offp, j);
        const R_xlen_t a1 = off_at(offp, j + 1);
        int b = pp[j];
        const int b1 = pp[j + 1];
        double s = 0.0;
        while (a < a1 && b < b1) {
            const int w_nr = wp[a];
            const int w_ob = ip[b];
            if (w_ob < w_nr) {
                ++b;
            } else if (w_ob > w_nr) {
                ++a;
            } else {
                s += xp[b];
                ++a;
                ++b;
            }
        }
        out[j] = s;
    }
    return sums;
}
