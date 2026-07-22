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
  vignette's "A Note On Computational Efficiency" section.

### Second pass (audit of the parallel cores themselves)

A profile of the first-pass cores exposed three structural costs, all
resolved within 0.12.0:

- **Bootstrap sampler complexity.** The conditional-binomial method is
  O(k_j) per replicate, but on wide corpora k_j ≫ N_j (measured: 2,500–6,500
  bins vs ~1,050 tokens per document). The sampler now chooses per document,
  from the data alone (thread-invariant): a Walker alias table draws the
  N_j tokens directly at O(N_j) per replicate when k_j > N_j, with the
  untouched bins contributing their closed-form Pearson mass; otherwise the
  binomial path runs with an early exit and suffix closed form once the
  remaining count reaches zero.
- **Serial share of the statistic pass.** At 4 threads the first-pass
  statistic was ~50% serial (Amdahl on 8.70 s → 3.47 s). The sparse-block
  densification now runs in parallel under the same slot pattern.
- **Sort length.** The envelope needs only the top-P_j head, so the full
  O(W log W) descending sort was replaced by `nth_element` selection with
  incremental slice sorts, O(W · widenings + P_j log P_j). Tie resolution is
  pinned to (probability descending, index ascending) — R's stable
  `order()` — closing the "partial sort" option recorded at 0.11.0 and
  making results platform-deterministic under tied fitted probabilities;
  the previously toolchain-sensitive Gibbs grid test now passes exactly
  everywhere and a dedicated tie test guards the contract.
- `optimal_topic()` additionally keeps the adapter extractions for the
  evaluation loop under a 256 MB grid budget instead of extracting every
  model twice.

### Third pass (the discrepancy-index family moves to compiled code)

The first audit had assessed the index family as "vectorized R whose cost
is dominated by BLAS" — measured at the scales the efficiency study now
targets, that no longer held. At medium-2 (J = 1,000, W = 10,000, grid of
10) the three-metric table took ~15 s and the partition ~3.5 s, dominated
by memory traffic over J × W dense temporaries (6–8 per metric per model,
plus logical→double coercions of the rare mask), not by the algebra;
extrapolated to the large scale this meant ~10 minutes and multi-GB
temporaries (`i_min` alone is 1.6 GB). Resolved within 0.12.0:

- **Fused index engine** (`src/index_core.cpp`). One compiled call per
  (model, block) computes every requested metric — SE, Pearson χ²,
  deviance — for the model and the no-topics baseline side simultaneously.
  At the document level the sparse counts are consumed directly: a dense
  pass accumulates the zero-count contributions and the collapsed-bin sums,
  a sparse pass over the nonzeros applies the corrections; the rare mask is
  bit-packed once per call (no logical→double coercion), and the baseline
  B = L_j·π_w is formed on the fly, never materialized. At the word level
  the fitted block from R BLAS is kept and the engine replaces the dense
  N/B temporaries. `optop_index_table()` therefore performs a single sweep
  per model for the whole metric set, and the D_null hoist (one baseline
  sweep per grid) is preserved. The support decomposition and the eps
  conventions are reproduced verbatim (χ² floors E element-wise before use
  and takes the min bin from floored sums; deviance keeps the 0·log 0
  convention and the unfloored E_min; SE is unfloored throughout) — the
  naive-reference tests remain the arbiter, and a golden check against the
  pre-engine implementation agreed to ≤ 4e-13 (document level) and
  ≤ 4e-17 (word level) on the fixtures.
- **Partition kernel** (`src/partition_core.cpp`).
  `optop_make_partition()` keeps the running minimum across the grid one
  vocabulary block at a time (one BLAS product per model per block feeding
  an in-place minimum) and writes the rare mask directly from the block
  buffer by exact comparison. The full J × W minimum matrix is never
  materialized; the mask is bit-identical to the R `pmin()` chain.
- **Envelope marshalling.** `optimal_topic()` now hands the statistic
  core's flattened envelope export (`bin_probs`, `bin_counts`) directly to
  the bootstrap core; the per-document list split survives only in the
  moment-matching path. This closes the "deferred, micro" transient noted
  in the second-pass benchmark (~1 GB on the large corpus).
- **API.** `n_threads` (default 1) on the three index functions,
  `optop_index_table()` and `optop_make_partition()`, with the established
  OpenMP contract: per-document/word private slots, fixed-order serial
  reduction, bit-identical output for any thread count (guarded in
  `test-parallel.R`).

**Methodological note (author's call, no code change).** With
τ_j = c / L_j, the harmonized support degenerates on wide, flat corpora:
as W grows at fixed document length, `min_K i_jw < c/L_j` eventually holds
for *every* word, C_j* covers the whole vocabulary, and the evaluation
support collapses to the single min bin. For such documents both D_K and
D_null reduce to `(N_min − E_min)²` with N_min = L_j and E_min ≈ L_j,
i.e. two O(ε²) rounding residuals, and the per-document
R² = 1 − D_K/D_null is a ratio of floating-point noise. Measured on the
benchmark corpora at the default c = 5: at the *small* scale
(J = 200, W = 2,000) 132 of 200 documents are fully collapsed but the
remaining 61 with D_null > 1 dominate the sums, so the micro indices stay
meaningful while the macro (equal-weight) average is already contaminated;
from medium-1 upward (W ≥ 5,000, L_j ≈ 600–1,500) **every** document
collapses — at medium-1, Σ_j D_null(SE) = 1.5e-20 — and the micro indices
themselves become noise ratios whose values legitimately differ between
algebraically equivalent implementations. The efficiency study therefore
evaluates the indices at c = 0.2 (no fully collapsed document at any
scale, median rare share 0.84–0.97), and its correctness gate compares
micro columns only. For the paper: if the indices are to be reported on
corpora in this regime, a guard may be warranted — a minimum support
size, a τ that adapts to W, or excluding (near-)collapsed documents from
the macro average.

## Deferred

- **Count-based statistic option**: the orthodox Pearson on counts
  (O = N_j·d, E = N_j·I), for which χ²_{P_j} is the textbook asymptotic —
  a different statistic than the published Eq. (8), so an alternative
  alongside it, not a replacement.

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
  makes this cleaner than it was at 0.9.7. *(Implemented in the 0.12.0
  second pass, see above.)*

### 0.12.0 parallel cores

Simulation study of `data-raw/benchmark-efficiency.R` (Intel Xeon @
2.10 GHz, 4 cores, Linux, gcc 13, reference BLAS, R 4.3.3; grids of 10 VEM
models, B = 200, q = 0.95; wall time from `system.time()` in a fresh R
process per configuration, peak RSS from GNU time). Thread sweeps verified
bit-identical; C++ vs R bootstrap p-values verified within Monte-Carlo
resolution.

Bootstrap-calibrated pass (statistic + calibration, whole grid):

| Corpus | R bootstrap (0.11.0) | C++ 1 thread | C++ 4 threads | total speedup |
|---|---|---|---|---|
| J = 200, W = 2,000 | 36.8 s | 26.1 s | 6.8 s | 5.4× |
| J = 500, W = 5,000 | 180.5 s | 120.6 s | 31.4 s | 5.7× |
| J = 1,000, W = 10,000 | 638.0 s | 426.5 s | 107.5 s | 5.9× |

Statistic-only pass:

| Corpus | 1 thread | 4 threads | speedup |
|---|---|---|---|
| J = 200, W = 2,000 | 0.34 s | 0.16 s | 2.2× |
| J = 500, W = 5,000 | 2.12 s | 0.90 s | 2.3× |
| J = 1,000, W = 10,000 | 8.70 s | 3.47 s | 2.5× |

Findings (first pass): the fused sampling+reduction is worth 1.4–1.5×
before threading; the bootstrap's document loop scales essentially linearly
(3.8–4.0× on 4 cores), the statistic sub-linearly (the gemm and the
sparse-block densification stay serial). Peak RSS is dominated by the
fitted models, is flat across implementations at the medium/large scales
(~12% saving at the small one), and does not grow with `n_threads`.

Second pass (adaptive sampler, envelope selection, parallel densification),
final study on four corpus scales — the former medium/large renamed
medium-1/medium-2, plus a J = 10,000 scale — with the thread sweep
{1, 2, 4} on the same 4-logical-core machine:

Bootstrap-calibrated pass (statistic + calibration, whole grid):

| Corpus | Models | R bootstrap (0.11.0) | C++ 1 thread | C++ 4 threads | total speedup |
|---|---|---|---|---|---|
| small: J = 200, W = 2,000 | 10 | 33.5 s | 18.1 s | 4.7 s | 7.1× |
| medium-1: J = 500, W = 5,000 | 10 | 165.9 s | 31.3 s | 8.7 s | 19.0× |
| medium-2: J = 1,000, W = 10,000 | 10 | 579.2 s | 67.0 s | 20.5 s | 28.3× |
| large: J = 10,000, W = 20,000 | 5 | 6,407 s (1 h 47 m) | 421.0 s | 141.9 s | 45.2× |

Statistic-only pass (1 thread → 4 threads): 0.37 → 0.18 s, 2.29 → 0.81 s,
9.36 → 3.17 s, 109.9 → 37.2 s (2.1–3.0×).

Findings:

- The pre-threading bootstrap gain grows monotonically with the vocabulary
  (1.9× / 5.3× / 8.6× / 15.2×), as expected from the O(N_j) alias path
  replacing the O(k_j) binomial path precisely where k_j ≫ N_j.
- Thread settings above the physical core count were measured once (6 and
  8 threads on this 4-core host) and sat within noise of the 4-thread
  times on every scale — pure oversubscription. They are omitted from the
  reported results, and the benchmark driver now caps the sweep at
  `parallel::detectCores()`, so wider grids appear only on machines that
  physically have the cores.
- Peak memory is governed by what must be held (the fitted models and, on
  the largest corpus, the exported envelope bins, ~7.5–8.6 GB): the two
  implementations stay within ~15% of each other in either direction, and
  memory does not grow with `n_threads`. On the largest corpus the
  compiled path sits ~14% above the R baseline because `.optop_boot_null()`
  flattens the envelope list before the core call — one transient copy;
  passing the core's flattened export directly would remove it. *(Done in
  the third pass, see above.)*

### 0.12.0 index engine (third pass)

Same harness, `indices` mode: `optop_make_partition()` plus a
three-metric, document-level `optop_index_table()` over the grid, on
counts reconstructed from the cached fits, at c = 0.2 (see the τ
degeneracy note for why the default c = 5 is unusable on these corpora;
the arithmetic volume does not depend on c). Baseline `indices_r` is the
pre-engine vectorized R implementation reimplemented verbatim in the
worker. Thread sweeps bit-identical; engine vs R indices agree to ≤ 4e-15
(micro and macro columns alike) at every scale.

| Corpus | Models | R implementation | engine 1 thread | engine 4 threads | total speedup | peak RSS |
|---|---|---|---|---|---|---|
| small: J = 200, W = 2,000 | 10 | 1.0 s | 0.18 s | 0.14 s | 7.1× | 345 → 274 MB |
| medium-1: J = 500, W = 5,000 | 10 | 6.4 s | 0.87 s | 0.65 s | 9.8× | 610 → 401 MB |
| medium-2: J = 1,000, W = 10,000 | 10 | 33.1 s | 3.7 s | 3.0 s | 11.1× | 1,216 → 625 MB |
| large: J = 10,000, W = 20,000 | 5 | 331.3 s | 39.8 s | 32.9 s | 10.1× | 9,260 → 3,269 MB |

Findings:

- The pre-threading gain (5.5× / 7.4× / 9.1× / 8.3×) is the fused single
  traversal itself: the R implementation's cost was the chain of J × W
  temporaries, which the engine never allocates. Within the total, the
  table sweep gains the most (large: 300 s → 21 s, 14.5×); the partition
  (large: 32 s → 12 s) remains BLAS-dominated.
- Thread scaling is modest by design of the workload (large: 1.2× from
  1 → 4 threads): the dense per-document pass is memory-bandwidth-bound
  and the partition's gemm stays serial inside the BLAS. `n_threads`
  helps, but the headline win is thread-free.
- Peak memory at the large scale drops 2.8× (9.3 → 3.3 GB): no `i_min`
  (1.6 GB), no dense count/baseline/difference blocks; what remains is
  the logical `rare_mask` (800 MB), the sparse counts, and the fits.

# Audit — GoF indices vs the revised theory paper (v0.13.0)

Formula-level audit of `optop_make_partition()`, `optop_make_baseline()`,
the `optop_index_*` family and `optop_index_table()` against Lewis &
Grossetti (2026), "Goodness-of-Fit Indices and Diagnostics for Topic
Models" (working paper, July 2026 revision).

## Verified correct (no change)

- Baseline `pi_glob` and `B_{j,min}` construction.
- Document-level SE, Pearson and deviance on the harmonized support
  `{non-rare} U {min}`, including the 0 log 0 convention and the unfloored
  min-bin expectation for the deviance.
- Macro = unweighted mean of per-document fits over `J+`.
- Word level on the unbinned vocabulary; word SE and Pearson; word Macro
  over `V+`; word Micro weights proportional to the baseline discrepancy.
- `tau_j = c / L_j`.

## Deviations found and fixed in 0.13.0

1. The harmonized union omitted the no-topics baseline; the running
   minimum of the partition kernel is now seeded at `pi_glob`.
2. The word-level fitted deviance lacked the Poisson linear correction
   `-(N - E)`; without it the word index could exceed one.
3. The Micro sums were not restricted to `J+`: degenerate documents
   leaked fitted discrepancy into the numerator. Degenerate units now
   carry `NA`.
4. The Pearson min-bin inclusion rule
   (`min(min_K E_min, B_min) >= c`, decided grid-wide, applied to both
   sides, reported) was absent; it now lives in the partition
   (`chisq_min_ok`) and the document kernel.
5. Features without paper backing deprecated: in-sample `ztest`
   (inference is reserved for held-out evaluation), `reopt` blending and
   `add_baseline_topic` (non-negativity enforcement contradicts the
   interpretation of negative indices as informative).
6. Default `c` moved from 5 to 1, the paper's adopted default for
   deviance-primary diagnostics.

Efficiency spot-run (medium-2 scale of the 0.12.0 study: J = 1,000,
W = 10,000, 10 synthetic models, 1 thread): partition 0.94 s, of which the
new rare-mass pass for the inclusion rule is 0.67 s; three-metric document
table 1.43 s; partition + table 2.37 s versus 3.7 s measured for the same
scale in the 0.12.0 study. The v0.12.0 gains are preserved.

## Roadmap (pure additions, each on its own branch)

- **Held-out machinery (paper §3.6)**: train/eval split, fold-in of the
  document-topic weights via `topicmodels::posterior(newdata)`,
  training-based harmonized partition and baseline-support conventions
  (out-of-support tokens to the min-bin, optional smoothing), per-document
  held-out fit, Macro confidence intervals, paired adjacent gains, the
  epsilon-adequacy selection rule, and the stabilized index variant.
- **Moment-based specification tests (paper §4)**: held-out residual
  moments against training-built instruments — the frequency-contrast
  screen (t-test), the frequency-strata Wald test, and the fit-stratified
  Wald test — with the sample-covariance Wald machinery and optional
  BH-adjusted marginal t-tests.
- **Cross-fitting helpers (paper Appendix B)** and standard errors for
  the Micro index and the Micro-Macro gap (delta method).

# Audit — WP4 revision (v0.14.0)

Formula-level audit of the GoF surface against the July 2026 revision of
Lewis & Grossetti, "Goodness-of-Fit Indices and Diagnostics for Topic
Models" (WP4).

## Verdict on the in-sample surface (0.13.0): conforms exactly

- Harmonized partition: WP4 restates C*_j as a union over the augmented
  grid K0 = K U {null}; identical to the baseline-seeded running minimum
  implemented in 0.13.0.
- WP4 now states explicitly that the min-bin exclusion rule is restricted
  to Pearson-type summaries because exhaustiveness is essential for the
  deviance KL bound: the kernel's chisq-only exclusion is confirmed.
- Deviance/Pearson/SE discrepancies, J+/V+ aggregation, Poisson word
  deviance, Micro weights, c = 1 default: all match.
- Notation only: fitted probabilities renamed i -> p across the docs and
  kernel comments.

## Implemented in 0.14.0 (the remaining theory)

- Held-out inference (paper 3.7): fold-in adapters, evaluation-side
  harmonized partition with the training baseline, out-of-support and
  smoothing conventions, per-document scores with the stabilized variant,
  Macro CI (Prop 2 i-iii), Micro and Micro-Macro gap delta-method SEs,
  paired adjacent gains and the epsilon-adequacy rule (Prop 2 iv, Def 1).
- Moment tests (paper 4): mean-zero training instruments (contrast,
  frequency strata, fit strata), Wald/t statistics with sample covariance
  (Prop 3), marginal t tests with adjustment.
- text2vec WarpLDA adapter via the optop_warplda() wrapper.

Defect found by the held-out tests and fixed: the doc kernel's
full-minus-rare accumulation cancelled catastrophically when a rare cell
carried a large count on an eps-floored expectation; non-rare cells now
accumulate directly (same cost, no cancellation).

## Roadmap

- Cross-fitting helpers (paper Appendix B): V-fold pooling with fold-aware
  variance estimators.
- The simulation harness of paper Section 5 lives outside the package by
  design (robustness testing of the theory, not package functionality).
