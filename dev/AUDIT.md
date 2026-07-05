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
bootstrap remains out of scope). The per-document bootstrap loop originally
lived in R on C-level primitives; since 0.12.0 it runs in
`src/calibration_core.cpp` under OpenMP (see below).

## Generalization to arbitrary engines (implemented in 0.11.0)

The 0.10.x core hard-coded `@gamma`/`@beta` S4 access and `optimal_topic()`
gated inputs with `is.LDA_VEM()` (rejecting even `LDA_Gibbs`). Resolved on
the `generalize-core` branch:

- **The C++ core takes plain matrices.** `optimal_topic_core(theta, phi,
  dfm_t, q, doc_map, return_envelope)` — no S4 slot access, no `exp()`, and
  the weighted dfm is transposed once in R for the whole grid instead of once
  per call. Adding a model class never touches C++ again: it is one R
  adapter method. As a by-product the hot path contains no R API calls,
  which is the thread-safety precondition for the OpenMP work below.
- **Adapter contract extended** to `list(theta, phi, K, docs, terms)` with a
  validator (`.optop_validate_theta_phi()`): simplex rows, unique mandatory
  document identifiers and terms, consistent dimensions. Methods:
  `TopicModel` (one workhorse covering topicmodels LDA VEM/Gibbs *and* CTM —
  `posterior()` is defined on the parent class and returns `exp(@beta)`
  verbatim), `textmodel_lda` (covers all three seededlda constructors, which
  share the class), `nlp_topic_fit` (NLPstudio: `dtw`/`tww`/`doc_ids`/`vocab`,
  with recursion into `model_object` when the weights are not stored), plus
  the informative `default` and the `LDA_t2v` stub.
- **Alignment hardened.** Two latent defects fixed: vocabulary alignment was
  never validated (a feature-permuted dfm was silently mis-scored — now
  reordered when it is the same set, error otherwise), and `doc_map` was
  built against model 1 and reused for the grid (now per model, with the
  document set intersected across all models). Alignment is identifier-based
  everywhere; positional alignment is never assumed.
- **Legacy frozen.** The pre-0.9.9 arithmetic lives verbatim in
  `src/optimal_topic_core_legacy.cpp`, used only by `selection = "legacy"`
  (still VEM-only), so bit-compatibility holds by construction; the file
  dies with the rule before v1.0.0.

## Parallelization (implemented in 0.12.0)

The OpenMP item deferred since 0.9.7 is closed, together with the bootstrap
port it grew into, on the `parallel-calibration` branch:

- **Statistic core.** The per-document sort/cumsum/Pearson loop — ~95% of
  the runtime since the blocked-gemm rewrite — runs under
  `#pragma omp parallel for` within each 256-document block. Each document
  writes a private result slot and the slots are reduced serially in
  document order, so the output is bit-identical for any `n_threads` and
  `Rcpp::checkUserInterrupt()` stays at the block level, outside the
  parallel region.
- **Bootstrap core** (`src/calibration_core.cpp`). The multinomial sampling
  (conditional binomial per bin) is fused with the Pearson reduction — the
  `bins × B` count and deviance matrices of the R implementation are never
  materialized — and parallelized over documents with the same
  slot-plus-serial-reduction pattern. R's RNG is not thread-safe, so the
  core owns its generator: each document draws from a `std::mt19937_64`
  seeded deterministically from `(seed, document index)` via splitmix64.
  Results are therefore bit-identical for any thread count and reproducible
  under `seed` (drawn once from the R session RNG when `NULL`, so
  `set.seed()` semantics are preserved). The draws are a different, equally
  valid stream than `rmultinom()`: calibrated p-values agree with 0.11.0 up
  to Monte-Carlo noise (`O(1/√B)`), not bit for bit — asserted against the
  pure-R oracle in `tests/testthat/test-parallel.R`.
- **API.** One `n_threads` argument on `optimal_topic()` (default 1),
  forwarded to both cores; `src/Makevars(.win)` adds
  `$(SHLIB_OPENMP_CXXFLAGS)` with an `#ifdef _OPENMP` fallback, so
  no-OpenMP toolchains (Apple clang) build and run single-threaded — those
  users keep the thread-free C++ speedup of the bootstrap, not the scaling.
- Benchmarks: see "0.12.0 parallel cores" in the benchmark section;
  reproducible via `data-raw/benchmark-efficiency.R`, reported in the
  vignette's "Computational efficiency" section.

## Deferred

- **Count-based statistic option**: the orthodox Pearson on counts
  (O = N_j·d, E = N_j·I), for which χ²_{P_j} is the textbook asymptotic —
  a different statistic than the published Eq. (8), so an alternative
  alongside it, not a replacement.

- **text2vec (WarpLDA) adapter**: the `LDA_t2v` stub still fails
  informatively; mapping `doc_topic_distr`/`topic_word_distribution` to the
  contract needs a verified normalization convention.

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

### 0.9.7 blocked-gemm rewrite (historical)

Timings of `optimal_topic(..., do_plot = FALSE)` on an Apple Silicon Mac
(R 4.6.1, Apple clang 21, `-O2`, reference BLAS), median of 3 runs, before
(commit `d44317f`, pre-rewrite) vs after the blocked-gemm rewrite. Outputs
agree to ~1e-14 in both settings.

| Corpus | Models | Before | After | Speedup |
|---|---|---|---|---|
| J = 1000, W = 4000 (VEM fits) | K ∈ {3, 5, 8, 12} | 1.41 s | 0.98 s | 1.4× |
| J = 2000, W = 8000 (synthetic weights) | K ∈ {5, 10, 20, 40, 80} | 6.43 s | 3.83 s | 1.7× |

At 0.9.7 the full grid's `gamma %*% exp(beta)` products accounted for only
~0.14 s of the 3.83 s in the large setting: after the rewrite, the
per-document *sorting* machinery (`sort_index` over W elements,
J × n_models times) dominated the runtime, not the algebra.

### 0.11.0 matrix core

Same settings re-run against the model-agnostic core (Apple Silicon Mac,
R 4.6.1, Apple clang 21, `-O2`, OpenBLAS; fresh synthetic corpora — the
0.9.7 harness was not preserved, so timings are indicative across audits,
exact within this one):

| Corpus | Models | Full pipeline | Bare core loop | gemm share |
|---|---|---|---|---|
| J = 1000, W = 4000 (VEM fits) | K ∈ {3, 5, 8, 12} | 0.83 s | — | — |
| J = 2000, W = 8000 (synthetic weights) | K ∈ {5, 10, 20, 40, 80} | 5.06 s | 5.11 s | 0.26 s |

Findings:

- **The generalization layer is free.** The full `optimal_topic()` pipeline
  (adapter pre-pass, contract validation, per-model document maps, feature
  checks) times the same as a bare loop over `optimal_topic_core()` calls on
  pre-adapted matrices — the R plumbing is not measurable next to the core.
- **The sort still dominates (~95%).** The gemm is ~5% of the large-setting
  runtime; the per-document `sort_index`/`cumsum` machinery is where the
  time goes, unchanged from the 0.9.7 profile. OpenMP over the document
  loop was therefore the highest-leverage follow-up — implemented in
  0.12.0, see the next benchmark.
- A further algorithmic option: the Eq. (8) statistic only needs the *set*
  of top-`P_j` words (sums are order-independent), so a
  partial-sort/selection scheme could replace the full descending sort;
  with the legacy pipeline frozen in its own translation unit, the modern
  core no longer needs to reproduce the `round(cum * 1e4)` cutoff, which
  makes this cleaner than it was at 0.9.7.
