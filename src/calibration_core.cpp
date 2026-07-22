#include <Rcpp.h>
#include <climits>
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
// Sampling strategy, chosen per document from the data alone (so the choice
// never depends on the thread count or the schedule):
// * k_j <= N_j: conditional binomial over the bins, O(k_j) per replicate,
//   with an early exit once the remaining count reaches zero — every
//   remaining bin then contributes its closed-form Pearson term p_b, held
//   in a precomputed suffix sum;
// * k_j > N_j (wide envelopes, the common case on large vocabularies): the
//   N_j tokens are drawn directly through a Walker alias table (O(k_j)
//   setup per document, O(N_j) per replicate). Counts are accumulated
//   sparsely; bins left untouched contribute sum(p_b) = 1 - sum over the
//   touched bins, so each replicate costs O(N_j) regardless of k_j.
//
// Threading and reproducibility contract:
// * each document owns a private RNG deterministically seeded from
//   (seed, j) via splitmix64, so its draws depend only on the seed and the
//   document index, never on the number of threads;
// * documents are processed in fixed-size blocks: the parallel region fills
//   one length-n_boot slot per document, and the slots are reduced serially
//   in document order, so the floating-point summation order is fixed and
//   the result is bit-identical for any n_threads;
// * no R API is touched inside the parallel region (R's RNG and allocator
//   are not thread-safe); interrupts are checked between blocks. The draws
//   are a different (equally valid) stream than R's rmultinom().

namespace {

// SplitMix64: the standard 64-bit seed scrambler (Steele, Lea & Flood 2014).
inline std::uint64_t splitmix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Walker alias table for O(1) draws from a discrete distribution.
struct AliasTable {
    std::vector<double> cut;
    std::vector<int> alias;

    void build(const double* p, const int k) {
        cut.assign(k, 0.0);
        alias.assign(k, 0);
        std::vector<double> scaled(k);
        std::vector<int> small, large;
        small.reserve(k);
        large.reserve(k);
        for (int b = 0; b < k; ++b) {
            scaled[b] = p[b] * k;
            if (scaled[b] < 1.0) small.push_back(b); else large.push_back(b);
        }
        while (!small.empty() && !large.empty()) {
            const int s = small.back(); small.pop_back();
            const int l = large.back(); large.pop_back();
            cut[s] = scaled[s];
            alias[s] = l;
            scaled[l] = (scaled[l] + scaled[s]) - 1.0;
            if (scaled[l] < 1.0) small.push_back(l); else large.push_back(l);
        }
        while (!large.empty()) { cut[large.back()] = 1.0; large.pop_back(); }
        while (!small.empty()) { cut[small.back()] = 1.0; small.pop_back(); }
    }

    inline int draw(std::mt19937_64& rng, const int k) const {
        const double u = std::generate_canonical<double, 53>(rng) * k;
        int b = static_cast<int>(u);
        if (b >= k) b = k - 1;
        const double frac = u - static_cast<double>(b);
        return (frac < cut[b]) ? b : alias[b];
    }
};

}  // namespace

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector optop_boot_null_core(const Rcpp::NumericVector& bin_probs,
                                         const Rcpp::IntegerVector& bin_counts,
                                         const Rcpp::NumericVector& doc_lengths,
                                         int n_boot,
                                         double seed,
                                         int n_threads,
                                         double doc_offset = 0)
{
    // the only R_xlen_t -> int narrowing in the package not backed by a
    // CSC or vocabulary bound: make the per-shard document cap explicit
    // rather than silent (a dgCMatrix shard cannot exceed it anyway)
    if (bin_counts.size() > static_cast<R_xlen_t>(INT_MAX)) {
        Rcpp::stop("the calibration shard holds more than 2^31 - 1 documents; "
                   "evaluate through optop_corpus() shards");
    }
    const int n_docs = static_cast<int>(bin_counts.size());
    if (doc_lengths.size() != n_docs) {
        Rcpp::stop("doc_lengths must have one entry per document");
    }
    if (n_boot < 1) {
        Rcpp::stop("n_boot must be positive");
    }
    if (doc_offset < 0) {
        Rcpp::stop("doc_offset must be nonnegative");
    }
    if (n_threads < 1) {
        n_threads = 1;
    }
    // sharded evaluation: the RNG stream of a document is keyed by its
    // GLOBAL index (shard offset + local index), so every document draws
    // bit-identical replicates whatever the sharding. The corpus replicate
    // T*_b is a sum of per-document terms; summing per-shard subtotals
    // groups the floating-point additions differently from the unsharded
    // serial reduction, so sharded totals agree to summation order (a few
    // ulp), not bit for bit.
    const std::uint64_t off64 = static_cast<std::uint64_t>(
        static_cast<long long>(doc_offset));

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
#pragma omp parallel num_threads(n_threads)
#endif
        {
            // per-thread scratch, sized on demand
            AliasTable at;
            std::vector<double> suffix;
            std::vector<int> count_scratch;
            std::vector<int> touched;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (int local = 0; local < block_len; ++local) {
                const int j = block_start + local;
                const double* p = probs + offset[j];
                const int k = counts[j];
                const long long n = static_cast<long long>(lengths[j] + 0.5);
                const double n_dbl = static_cast<double>(n);
                const double k_dbl = static_cast<double>(k);
                double* T_doc =
                    contrib_ptr + static_cast<std::size_t>(local) * n_boot;

                // private, thread-count-independent stream for this document,
                // keyed by the global document index
                std::mt19937_64 rng(
                    splitmix64(seed64 + off64
                               + static_cast<std::uint64_t>(j) + 1));

                if (k > n) {
                    // wide envelope: draw the n tokens directly, accumulate
                    // counts sparsely, and add the closed-form contribution
                    // of the untouched bins
                    at.build(p, k);
                    double p_total = 0.0;
                    for (int bin = 0; bin < k; ++bin) {
                        p_total += p[bin];
                    }
                    if (static_cast<int>(count_scratch.size()) < k) {
                        count_scratch.assign(k, 0);
                    }
                    touched.clear();
                    touched.reserve(static_cast<std::size_t>(n));

                    for (int b = 0; b < n_boot; ++b) {
                        for (long long t = 0; t < n; ++t) {
                            const int bin = at.draw(rng, k);
                            if (count_scratch[bin]++ == 0) {
                                touched.push_back(bin);
                            }
                        }
                        double pearson = 0.0;
                        double touched_mass = 0.0;
                        for (const int bin : touched) {
                            const double pb = p[bin];
                            const double dev =
                                static_cast<double>(count_scratch[bin]) / n_dbl - pb;
                            pearson += dev * dev / pb;
                            touched_mass += pb;
                            count_scratch[bin] = 0;
                        }
                        touched.clear();
                        // every untouched bin contributes (0 - p_b)^2 / p_b
                        pearson += p_total - touched_mass;
                        T_doc[b] = k_dbl * pearson;
                    }
                } else {
                    // narrow envelope: conditional binomial over the bins,
                    // with the suffix closed form once the count is exhausted
                    suffix.assign(k + 1, 0.0);
                    for (int b = k - 1; b >= 0; --b) {
                        suffix[b] = suffix[b + 1] + p[b];
                    }

                    for (int b = 0; b < n_boot; ++b) {
                        long long remaining = n;
                        double prob_left = 1.0;
                        double pearson = 0.0;

                        for (int bin = 0; bin < k; ++bin) {
                            if (remaining == 0) {
                                // all further bins draw zero counts: their
                                // Pearson terms sum to the suffix mass
                                pearson += suffix[bin];
                                break;
                            }
                            long long c;
                            if (bin == k - 1) {
                                // the last bin takes whatever is left: the
                                // bin probabilities sum to 1 per document
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
                            const double dev =
                                static_cast<double>(c) / n_dbl - p[bin];
                            pearson += dev * dev / p[bin];
                        }

                        T_doc[b] = k_dbl * pearson;
                    }
                }
            }
        }

        // serial reduction in document order: the summation order (and hence
        // the floating-point result) never depends on the schedule
        for (int local = 0; local < block_len; ++local) {
            const double* T_doc =
                contrib_ptr + static_cast<std::size_t>(local) * n_boot;
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
