# Audit — `optimal_topic()` / `src/optimal_topic_core.cpp`

Status of the findings from the performance and correctness audit of the
optimal-topic pipeline. This file tracks what was found, what has been fixed on
the `optimal-topic-core` branch, and what is deliberately deferred.

## Performance (fixed on this branch)

1. **Dense W×K temporary per document.** The per-document hot loop built
   `tww.each_row() % dtw.row(j)` — a dense W×K matrix — and row-summed it.
   That row sum is exactly row *j* of the matrix product
   `gamma %*% exp(beta)`, so the whole inner-loop algebra collapses to one
   BLAS `gemm` per (model, document-block). What remains per document is only
   a descending sort, a cumulative sum, and a scalar chi-square reduction.
2. **Element-wise reads from a sparse row.** Each document row of the weighted
   dfm was copied to a dense vector element by element from an
   `arma::sp_mat` (CSC storage, so row access is slow), and this dense copy
   was redone for *every model* instead of once. The rewrite densifies the
   dfm in document blocks, once per block, shared across all models.
3. **Blocked processing.** Both the dense dfm copy and the `gemm` output are
   materialized in fixed-size document blocks rather than as full J×W
   matrices, so memory stays bounded for large corpora.
4. Dead parameters `docs` and `n_features` were removed from the C++
   signature; the J×4 `regstats` matrix was replaced by two running
   accumulators (only the column sums were ever used).

## Correctness

### Fixed on this branch

- **Document-order alignment.** `optimal_topic()` checked only *membership*
  of `docid(weighted_dfm)` in `lda_models[[1]]@documents`, never *order*.
  Row *j* of the dfm was assumed to correspond to row *j* of `@gamma`; a dfm
  whose rows were permuted relative to the fitted models was silently
  mis-scored. Each dfm row is now paired with its `@gamma` row by document
  identifier (a `doc_map` index vector passed to the C++ core), and the core
  stops if the mapping points past the rows of `@gamma`.

### Resolved in 0.9.9 — calibrated to the published test (Test 1, Eq. 8)

Reading the JMLR paper (Lewis & Grossetti 2022, 23(58)) settled both flags and
uncovered three more deviations; all five are fixed on the `test-calibration`
branch, with the old pipeline frozen behind `selection = "legacy"`
(deprecated, bit-identical to v0.9.8, removal before v1.0.0):

- **Tail and df.** The paper claims `OpTop ~ χ²` with `df = Σ_j P_j` and
  frames adequacy as "unable to reject" — an upper-tail test on the *raw*
  statistic. The code used the lower tail with `df = 1` on the *standardized*
  statistic. Now: `pchisq(raw, df = Σ P_j, lower.tail = FALSE)`; the reported
  `OpTop` column stays standardized (`raw / df`, the paper's Figure 2 scale).
- **Selection.** The paper selects the global minimum of the standardized
  statistic and never uses p-values for selection; the old "first
  `pval ≤ alpha`" lower-tail rule was an implementation-era invention. Now:
  `selection = c("sequential", "min", "legacy")` — the default sequential
  adequacy scan (smallest K with `pval > alpha`, global-minimum fallback),
  the paper's `"min"`, and the frozen legacy rule.
- **Per-document scaling.** Eq. (8) scales each document's Pearson term by
  `(P_j + 1)` (the number of bins); the code used `P_j`.
- **Cutoff.** Footnote 5 collapses the least-probable words with total mass
  strictly below `I^K = 0.05`; the code rounded the cumulative mass to 4
  decimals, put the crossing word in the tail (collapsed mass *exceeding*
  `1 − q`), and defaulted to `q = 0.80`. Now: exact comparison, crossing word
  kept, default `q = 0.95 = 1 − I^K`.
- **Calibration caveat (documented, not "fixed").** With `df = Σ P_j` in the
  thousands, χ² p-values saturate near 0/1 unless the fit is genuinely
  borderline; when the sequential rule degenerates, the minimum of the
  standardized curve carries the information — which is why the paper's case
  study uses the min rule. Establishing a sharper null law (asymptotics or
  simulation calibration) remains future methodological work for the paper.

## Null calibration for Test 1 (implemented in 0.10.0)

The χ² reference of Eq. (8) is a yardstick, not an exact null law: the
classical Pearson asymptotics require the statistic to be scaled by the
document *length* `N_j` (counts), whereas Eq. (8) works on proportions scaled
by the bin count `P_j + 1` — a per-document scale gap of roughly
`N_j / (P_j + 1)` — and the expected probabilities are estimated by VEM from
the same data. Hence the saturated p-values documented in `?optimal_topic`.

Implemented via `optimal_topic(calibrate = ...)` (see the Calibration section
of the man page for the full reasoning):

1. **`"bootstrap"`** — parametric bootstrap under the conditional
   fitted-model null. Because the envelope is data-free and a multinomial
   collapsed over bins is multinomial on the collapsed probabilities, null
   replicates are drawn directly on the `P_j + 1` bins exported by the C++
   core (`return_envelope`), at ~`n_boot × df` flops per model; empirical
   p-value `(1 + #{T* ≥ T})/(n_boot + 1)`.
2. **`"moment"`** — exact Haldane (1937) multinomial moments of each
   document's Pearson term, summed and Satterthwaite-matched to `a·χ²_ν`.
   Closed form; validated in the test suite against the bootstrap's empirical
   moments.

Both are conditional on the fitted Θ̂, Φ̂ (no per-replicate refit — the double
bootstrap remains out of scope). The per-document bootstrap loop lives in R
on C-level primitives; it is the designated C++/OpenMP port if corpora with
J ≫ 10⁴ ever make it the bottleneck.

## Deferred

- **Count-based statistic option**: the orthodox Pearson on counts
  (O = N_j·d, E = N_j·I), for which χ²_{P_j} is the textbook asymptotic —
  a different statistic than the published Eq. (8), so an alternative
  alongside it, not a replacement.

- **OpenMP over the document loop.** Correct parallelization axis, but only
  worth adding after the `gemm` rewrite is benchmarked: it requires a
  `src/Makevars` with `$(SHLIB_OPENMP_CXXFLAGS)` and keeping
  `Rcpp::checkUserInterrupt()` out of the parallel region. See the benchmark
  section — revisit only if the per-document sort/cumsum still dominates.
- **Generalizability.** The core hard-codes `@gamma` / `@beta` S4 access and
  `optimal_topic()` gates inputs with `is.LDA_VEM()` (which even rejects
  `LDA_Gibbs`). The right abstraction already exists in `R/utils.R`
  (`optop_as_theta_phi()`): extract theta/phi in R and pass plain matrices to
  C++, with new adapter methods for `CTM_VEM`, seededlda's `textmodel_lda`,
  and NLPstudio's `nlp_topic_fit` (dtw → theta, tww → phi). Planned as a
  follow-up branch once the efficiency work is merged.

## Vignette build design (CRAN-ready)

The vignette's `optimal_topic()` calls evaluate **live** when the local model
cache (`data-raw/VEM_models_inaugural.rds`, ~170 MB, git- and build-ignored)
is present — the maintainer's `R CMD build` therefore ships an HTML with
genuine console output — and fall back to **replaying the captured console
streams** stored in the shipped 10 KB results bundle when the cache is absent
(CI checkouts, CRAN's re-building of vignette outputs inside the pruned
tarball). Both paths produce the same numbers because the bundle is derived
from the same cached fits. This conditional-evaluation pattern is standard,
CRAN-sanctioned practice for expensive vignettes; if CRAN's rebuild
environment ever becomes a problem, the fallback plan is the rOpenSci
"precompiled vignette" pattern (`OpTop.Rmd.orig` knitted locally into a fully
static `OpTop.Rmd`).

## Benchmark

Timings of `optimal_topic(..., do_plot = FALSE)` on an Apple Silicon Mac
(R 4.6.1, Apple clang 21, `-O2`, reference BLAS), median of 3 runs, before
(commit `d44317f`, pre-rewrite) vs after the blocked-gemm rewrite. Outputs
agree to ~1e-14 in both settings.

| Corpus | Models | Before | After | Speedup |
|---|---|---|---|---|
| J = 1000, W = 4000 (VEM fits) | K ∈ {3, 5, 8, 12} | 1.41 s | 0.98 s | 1.4× |
| J = 2000, W = 8000 (synthetic weights) | K ∈ {5, 10, 20, 40, 80} | 6.43 s | 3.83 s | 1.7× |

**Where the time goes now.** The full grid's `gamma %*% exp(beta)` products
account for only ~0.14 s of the 3.83 s in the large setting: after the
rewrite, the per-document *sorting* machinery (`sort_index` over W elements,
J × n_models times) dominates the runtime, not the algebra. Consequences:

- The headline gain of the rewrite is asymptotic in K (the removed temporary
  was W × K per document) and in sparsity handling; at small K the sort was
  already the bottleneck.
- **OpenMP over the document loop is now the highest-leverage follow-up**:
  documents are independent and the per-model reductions are two scalars, so
  the sort-dominated loop should scale near-linearly with cores.
- A further algorithmic option: the chi-square statistic only needs the
  *set* of top-`icut` words (sums are order-independent), so a
  partial-sort/selection scheme could replace the full descending sort —
  but the `round(cum * 1e4)` cutoff must be reproduced exactly, so this
  needs care.
