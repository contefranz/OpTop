#include <Rcpp.h>
#include <cstdint>
#include <random>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Parallel parametric bootstrap of the Test 1 null (see R/calibration.R for
// the statistical background). The conditional fitted-model null draws each
// document as Multinomial(N_j, p_j) directly on its collapsed envelope bins
// and accumulates T* = sum_j k_j * sum_b (c_jb/N_j - p_jb)^2 / p_jb, one
// value per replicate.
//
// Threading and reproducibility contract:
// * each document owns a private RNG deterministically seeded from
//   (seed, j) via splitmix64, so its draws depend only on the seed and the
//   document index — never on the number of threads or the schedule;
// * documents are processed in fixed-size blocks: the parallel region fills
//   one length-n_boot slot per document, and the slots are reduced serially
//   in document order, so the floating-point summation order is fixed and
//   the result is bit-identical for any n_threads;
// * no R API is touched inside the parallel region (R's RNG and allocator
//   are not thread-safe); interrupts are checked between blocks. The
//   multinomial is drawn by the conditional binomial method through
//   std::binomial_distribution — a different (equally valid) stream than
//   R's rmultinom().

// SplitMix64: the standard 64-bit seed scrambler (Steele, Lea & Flood 2014).
static inline std::uint64_t splitmix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector optop_boot_null_core(const Rcpp::NumericVector& bin_probs,
                                         const Rcpp::IntegerVector& bin_counts,
                                         const Rcpp::NumericVector& doc_lengths,
                                         int n_boot,
                                         double seed,
                                         int n_threads)
{
    const int n_docs = bin_counts.size();
    if (doc_lengths.size() != n_docs) {
        Rcpp::stop("doc_lengths must have one entry per document");
    }
    if (n_boot < 1) {
        Rcpp::stop("n_boot must be positive");
    }
    if (n_threads < 1) {
        n_threads = 1;
    }

    // per-document offsets into the flattened bin-probability vector
    std::vector<std::size_t> offset(n_docs + 1, 0);
    for (int j = 0; j < n_docs; ++j) {
        if (bin_counts[j] < 1) {
            Rcpp::stop("every document needs at least one envelope bin");
        }
        offset[j + 1] = offset[j] + static_cast<std::size_t>(bin_counts[j]);
    }
    if (offset[n_docs] != static_cast<std::size_t>(bin_probs.size())) {
        Rcpp::stop("bin_probs and bin_counts disagree on the envelope size");
    }

    const double* probs = &bin_probs[0];
    const int* counts = &bin_counts[0];
    const double* lengths = &doc_lengths[0];
    const std::uint64_t seed64 =
        static_cast<std::uint64_t>(static_cast<long long>(seed));

    const int block_size = 1024;
    std::vector<double> T_null(n_boot, 0.0);
    std::vector<double> contrib;  // block_len x n_boot, one row per document

    for (int block_start = 0; block_start < n_docs; block_start += block_size) {
        Rcpp::checkUserInterrupt();

        const int block_len = std::min(block_size, n_docs - block_start);
        contrib.assign(static_cast<std::size_t>(block_len) * n_boot, 0.0);
        double* contrib_ptr = contrib.data();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
        for (int local = 0; local < block_len; ++local) {
            const int j = block_start + local;
            const double* p = probs + offset[j];
            const int k = counts[j];
            const long long n = static_cast<long long>(lengths[j] + 0.5);
            const double n_dbl = static_cast<double>(n);
            const double k_dbl = static_cast<double>(k);
            double* T_doc = contrib_ptr + static_cast<std::size_t>(local) * n_boot;

            // private, thread-count-independent stream for this document
            std::mt19937_64 rng(
                splitmix64(seed64 + static_cast<std::uint64_t>(j) + 1));

            for (int b = 0; b < n_boot; ++b) {
                long long remaining = n;
                double prob_left = 1.0;
                double pearson = 0.0;

                for (int bin = 0; bin < k; ++bin) {
                    long long c;
                    if (bin == k - 1) {
                        // the last bin takes whatever is left: the bin
                        // probabilities sum to 1 per document
                        c = remaining;
                    } else {
                        double pr = p[bin] / prob_left;
                        if (pr > 1.0) pr = 1.0;
                        if (pr < 0.0) pr = 0.0;
                        std::binomial_distribution<long long> draw(remaining, pr);
                        c = draw(rng);
                        remaining -= c;
                        prob_left -= p[bin];
                        if (prob_left < 0.0) prob_left = 0.0;
                    }
                    const double dev = static_cast<double>(c) / n_dbl - p[bin];
                    pearson += dev * dev / p[bin];
                }

                T_doc[b] = k_dbl * pearson;
            }
        }

        // serial reduction in document order: the summation order (and hence
        // the floating-point result) never depends on the schedule
        for (int local = 0; local < block_len; ++local) {
            const double* T_doc = contrib_ptr + static_cast<std::size_t>(local) * n_boot;
            for (int b = 0; b < n_boot; ++b) {
                T_null[b] += T_doc[b];
            }
        }
    }

    return Rcpp::NumericVector(T_null.begin(), T_null.end());
}

//' @keywords internal
// [[Rcpp::export]]
bool optop_openmp_available()
{
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}
