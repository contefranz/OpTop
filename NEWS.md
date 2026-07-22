# OpTop 0.19.0

The cross-fitting release, closing the last statistical gap between the
package and the paper. The API surface is now frozen for 1.0.0.

### Renamed

* **`optimal_topic()` is now `optop_select()`.** The last function outside
  the package-wide `optop_` naming joins it before the 1.0.0 freeze, so
  typing `optop_` completes the entire API. The old name keeps working as
  a delegating alias that notes the rename once per session and will not
  be removed before a 2.0; scripts need no immediate change.

### New features

* **`optop_crossfit(dtm, K, fit_fun, V = 5, ...)`** implements the V-fold
  cross-fitted indices of the paper's Appendix B. The documents are dealt
  into length-balanced folds, the model grid is refit on each fold's
  complement by a user-supplied fitting function, and every fold is
  scored held-out through the existing machinery, so each document is
  evaluated exactly once by models that never saw it. Point estimates
  pool the out-of-fold discrepancies at the document level under the
  usual floor and stabilization conventions; standard errors for Macro,
  Micro, the Micro-Macro gap, and the adjacent gains are fold-replicate
  (cluster) estimates over the V folds, so the reported uncertainty
  includes the fold-to-fold variation induced by refitting. The
  epsilon-adequacy rule runs on the cross-fitted gains with these
  bounds. Any engine with an adapter and a fold-in path qualifies;
  sharded corpora are out of scope by design.
* `print()` and `plot()` methods for the cross-fitted result, sharing
  the holdout views: the pooled Macro curve with its fold-replicate
  confidence band and the cross-fitted adjacent-gains view.

### Documentation

* New vignette subsection **"V-Fold Cross-Fitting"** inside the held-out
  section, picking up exactly the two caveats of the single split (refit
  variability and split arbitrariness) and demonstrating the remedy on
  the inaugural corpus with a cached five-fold refit.

# OpTop 0.18.0

The interface release: everything a user touches at the console gets a
face. New model entry points, classed results with compact print methods,
a plot method for held-out runs, and the two vignette sections the
feature set had outgrown.

### New features

* **`optop_model(theta, phi)`**: the public entry point for engines the
  package has no adapter for. Any software that estimates a
  document-topic and a topic-word matrix now enters the toolchain through
  a validated constructor (named documents and terms, nonnegative
  weights, unit row sums) and flows through the in-sample tools like a
  native fit. Fold-in refuses it with a pointed message: no engine sits
  behind it, so held-out evaluation wants the engine's own object.
* **`optop_align_dtm_to_models()` is now exported.** The remediation
  helper the alignment errors have pointed at since 0.13.0 ships with a
  reference page, a runnable example, and direct tests.
* **Every result prints itself.** Partitions, the three discrepancy-index
  results, and the grid table now carry classes next to their unchanged
  list (or `data.table`) structure, joining the already classed holdout
  and moment-test results; compact cli print methods cover all of them
  plus the model container. Dollar-access code is unaffected.
* **`plot(holdout)`**: the Macro curve with its confidence band across
  metrics, and a `which = "gains"` view drawing the adjacent paired gains
  with one-sided upper bounds, the epsilon line, and the K selected by
  the adequacy rule.

### Documentation

* New vignette section **"The Null-Discrepancy Floor"**: why the Macro
  denominator needs a floor, the default `min_null = c`, and the real
  collapsed speech in the inaugural corpus whose single unfloored ratio
  drags the Macro index while Micro barely moves.
* New vignette section **"Other Engines and Bare Matrices"**: the adapter
  layer, `optop_warplda()` with a VEM-versus-WarpLDA comparison on one
  scale, `optop_model()`, and the alignment helper.

### Internal

* The dense-row reconstructor orphaned by the 0.16.0 teardown is
  removed; the `optop_fold_in()` methods gain direct per-engine contract
  tests.

# OpTop 0.17.0

The kernel release: the compiled cores lose their last avoidable
large-vocabulary costs, every long-running kernel becomes interruptible,
and the package website gains an operations article for very large
corpora. One numerical change, documented below.

### Performance

* **The collapsed-tail term of Test 1 no longer scans the vocabulary.**
  For row-stochastic model weights and row-normalized proportions, the two
  document totals in the tail term are identically 1, so the tail mass is
  now the exact complement of the head mass. This removes an O(J x W) pass
  per model, the dominant non-BLAS cost of `optimal_topic()` at large
  vocabularies. Together with the boundary change below, the scale
  benchmark's `optimal_topic()` rows run 2.2x to 2.6x faster than 0.15.0
  (one million documents by 50k terms: 1617 s down to 739 s; see
  `dev/benchmark-scale.md`).
* **Envelope search boundary.** The head selection starts at a fixed 1024
  candidates instead of a quarter of the vocabulary; the strict total
  order of the tie-breaking comparator guarantees an identical envelope,
  found with far less sorting at large W.
* **Block buffers are allocated once per call**: the fitted-probability
  block product (about 100 MB per 256-document block at W = 50k) now
  writes into preallocated memory through no-copy views.
* Doc-kernel gemv batching was evaluated and rejected: the merge-join
  kernel runs three orders of magnitude below the pipeline bottleneck and
  a gemv would only add a gather copy (`dev/benchmark-scale.md`).

### Breaking changes

* **One numerical change**: replacing the summed tail totals with their
  exact complements shifts results by at most about 1e-12 relative on
  adversarial fixtures (the exact constant replaces a W-term floating
  sum, so the new convention is the more accurate one). All reference
  tolerances hold unchanged; calibrated p-values can move by one
  bootstrap resolution step when a replicate sits within ulps of the
  observed statistic. The 0.15.0 Pearson min-bin note set the precedent
  for this kind of documented bit-level change.

### Robustness

* **Every partition kernel is now interruptible**: the six kernels iterate
  inside serial macro-chunks with an interrupt check between chunks, the
  pattern the index and calibration kernels already used, so multi-minute
  passes over tens of millions of documents can be cancelled from the
  console. Chunking regroups no arithmetic; results are bit-identical.
* The bootstrap core now errors explicitly past 2^31 - 1 documents per
  shard (previously an undiagnosed narrowing), pointing at
  `optop_corpus()`.

### New features

* `options(optop.cache_mb = )` configures the extracted-weights cache
  budget of the evaluation loops (default 256 MB, both the plain and the
  sharded path). Results never depend on the budget, only wall time.
* New package website article **"OpTop at very large scale"**: sharding
  with `optop_corpus()`, lazy model loaders, the cache budget, thread and
  BLAS configuration (including the `mclapply` oversubscription trap), and
  memory planning with reference timings. The article is website-only and
  never builds on CRAN.

# OpTop 0.16.0

The teardown release on the road to v1.0.0: every deprecated surface leaves
the package, three hard dependencies retire, and the compiled kernels gain
a sanitizer gate in continuous integration. No statistical behavior of the
surviving paths changes.

### Breaking changes

* **The pre-0.9.8 stability API is removed** after its long deprecation:
  `topic_stability()`, `agg_topic_stability()`, `agg_document_stability()`,
  `get_topic_models()`, and the internal `topic_match()` are gone, together
  with their compiled cores (`topic_match_core.cpp`, `normalize_columns`).
  The discrepancy indices and the held-out tools are their replacement.
* **`selection = "legacy"` is removed** from `optimal_topic()` (deprecated
  since 0.9.9). The argument now admits `"sequential"` and `"min"`; the
  frozen pre-0.9.9 translation unit is deleted, retiring the package's last
  `arma::sp_mat` crossing and the lower-tail p-value convention.
* **The 0.13.0-deprecated index arguments are removed**: `ztest`, `reopt`,
  and `add_baseline_topic` leave `optop_index_se()`, `optop_index_chisq()`,
  `optop_index_deviance()`, and `optop_index_table()`, together with the SE
  re-optimization path and the in-sample Z-test helper. The methodology
  reserves inference for held-out evaluation (`optop_index_holdout()`,
  `optop_moment_test()`).
* **`optimal_topic(lda_models = )` is removed** (renamed to `topic_models`
  in 0.11.0), and `optop_make_partition()` drops the `block` argument that
  had been inert since 0.15.0.
* **`sim_dfm()` is self-contained**: the generator now draws the LDA
  Dirichlet-multinomial corpus inside OpTop instead of delegating to
  LDATS. The simulated distribution is unchanged, but the draws for a
  given `seed` differ from earlier versions.

### Dependencies

* `LDATS`, `patchwork`, and `lifecycle` leave `Imports`, and `withr`
  leaves `Suggests`; each had lost its last caller with the removals
  above. No new dependency enters.

### Infrastructure

* New `sanitizers` workflow in continuous integration: the test suite runs
  under AddressSanitizer and UndefinedBehaviorSanitizer (r-hub clang-asan
  container) on kernel-touching changes and weekly, failing on any
  sanitizer finding.
* Build hygiene: stray `.Rhistory` files and rendered vignette artifacts
  are excluded from the build, and `sim_dfm()` gains a dedicated test file
  plus the package's first runnable example.

# OpTop 0.15.0

The scale release: the package now evaluates corpora of tens of millions of
documents on server machines with bounded memory, exactly. Three
architectural changes carry it, each gated by equivalence tests against the
previous behavior.

### New Features

* **Sharded corpora** (`optop_corpus()`): a dtm, a list of document shards,
  or per-shard file paths with a reader are accepted wherever a dtm or
  weighted dfm is accepted, by `optimal_topic()`, the partition and baseline
  builders, the index family and grid table, the held-out layer, and the
  moment tests. Every statistic is a sum of independent per-document terms,
  so sharded evaluation is exact: per-document outputs are bit-identical to
  an unsharded run and concatenate, sums agree up to floating-point
  summation order, and the calibration bootstrap keys each document's random
  stream by its global index, so a one-shard corpus reproduces the plain
  call bit for bit, calibrated p-values included. Sharding is how the
  package passes the container cap of R sparse matrices (2^31 - 1 nonzeros
  per dgCMatrix): shards stay under it individually while the corpus does
  not.

* **Lazy models**: elements of a model grid may be loader functions
  returning a fit; `optimal_topic()`, the partition builder, the index
  family, the holdout, and the moment tests materialize one model at a time
  and release it, so the grid never needs to sit in memory at once.
  Extracted weight matrices are still cached under a fixed budget when they
  fit.

### Performance

* **Sparse harmonized partition**: the partition stores the complement of
  the rare set, the per-document non-rare word lists, bounded in total by
  tokens divided by `c`. The construction is candidate-based and exact:
  because the null baseline belongs to the augmented union, a cell can be
  non-rare only if its baseline probability clears `tau_j`, so each
  document's candidates form a prefix of the baseline-sorted vocabulary.
  Cost falls from `O(J * W * sum_K K)` dense products to
  `O((tokens / c) * sum_K K)`, and no J x W object of any kind exists
  anymore (the dense `rare_mask`, its packed bits, and the per-document
  full-vocabulary traversal are gone). On 100,000 documents by 30,000
  features the partition builds in 0.2 s as a 53 MB object where the dense
  mask alone was 12 GB.

* **Merge-join document kernel**: each document is scored by merge-joining
  its non-rare list with its sparse counts at `O(|NR_j| * K + nnz_j)`; the
  word kernel blocks documents internally. The corpus crosses the C++
  boundary everywhere as zero-copy dgCMatrix slot views: no Armadillo
  sparse container, no per-model copy, no shape limit.

### Breaking Changes

* **Partition format 2**: `optop_make_partition()` returns the sparse
  fields (`nonrare_offsets`, `nonrare_words`, `vocab`, `format = 2L`)
  instead of the dense `rare_mask`. Partitions saved by 0.13.0 to 0.14.3
  are upgraded automatically on first use, with an alert; the `block`
  argument is accepted and ignored. The deprecated SE re-optimization path
  reconstructs its dense rows from the sparse lists and remains a
  small-corpus path.

* **Pearson min-bin convention**: the fitted min-bin mass is now the
  unfloored sum of the rare fitted probabilities, floored once at the
  numerical epsilon, unifying it with the null side (previously the fitted
  side summed element-wise floored expectations, which a complement-based
  construction cannot reproduce). The difference is bounded by `W * eps`
  in absolute terms on a bin that enters the Pearson family only when its
  mass clears `c`; reference implementations follow the same convention.

# OpTop 0.14.3

A large-corpus indexing release.

### Fixes

* **64-bit Armadillo indexing** (`-DARMA_64BIT_WORD` in both Makevars): the
  compiled cores receive the document-term matrix as an Armadillo sparse
  matrix, whose constructor requires the documents-by-features shape to be
  addressable by its integer index type. With the previous 32-bit indices any
  corpus whose `J * W` exceeded `2^31 - 1` cells failed at the boundary with
  `SpMat::init(): requested size is too large`, regardless of how sparse the
  matrix actually was (observed on 3.1 million documents by 57,440 features).
  All cores now compile with 64-bit indices, which also protects the dense
  document-topic matrices whose element count passes `2^31` at extreme scale.

### Known limits at very large scale

* A single `dgCMatrix` (and hence a single quanteda dfm) holds at most
  `2^31 - 1` nonzero entries, an R container limit no compiler flag lifts;
  corpora beyond it need document sharding, planned for the next minor
  release together with a zero-copy sparse boundary and a sparse partition
  representation. Until then `optimal_topic()` is the supported entry point
  at multi-million-document scale: `optop_make_partition()` and the index
  family still materialize dense documents-by-features objects internally
  and remain intended for corpora of moderate size.

# OpTop 0.14.2

A documentation release. The README now states the per-platform build
prerequisites for installing from source (Rtools on Windows; the Command Line
Tools and the official R gfortran toolchain on macOS, the latter required to
link the BLAS/LAPACK Fortran runtime; the compiler, gfortran and R headers on
Linux), and adds a note on cross-platform reproducibility: OpTop is
deterministic for a given set of fitted models, while the small differences
seen across operating systems originate in the upstream VEM fit. No code
changes.

# OpTop 0.14.1

The regularity release: a family-generic null-discrepancy floor protects the
document-level aggregations from structurally degenerate documents.

### Breaking Changes

* **Null-discrepancy floor** (`min_null`, default: the partition constant
  `c`): a document enters the Micro and Macro aggregations of
  `optop_index_se()`, `optop_index_chisq()`, `optop_index_deviance()`,
  `optop_index_table()` and `optop_index_holdout()` only when its baseline
  discrepancy satisfies `D_j(null) >= min_null`, uniformly across the three
  families. At partition resolution `tau_j = c/L_j` a document whose support
  collapses into the min bin has observed equal to baseline counts by
  construction, so `D_j(null) ~ 0` in every family at once and the
  document-level index `1 - D_j(K)/D_j(null)` is unbounded; held-out
  simulation evidence showed a handful of such documents destroying the
  Macro curves and their confidence intervals while the discrepancy-weighted
  Micro index was immune. The floor sharpens the conditioning event of the
  held-out inference from `{D_j(null) > 0}` to `{D_j(null) >= c}`; the
  estimand becomes the average fit among documents whose baseline
  discrepancy is at least the resolution of the evaluation support.
  Excluded documents carry an `NA` index, are counted per call
  (`n_null_excluded`, plus a share in the in-sample results and a summary
  column in the held-out table), and are announced once per metric.
  `min_null = 0` restores the previous strict-positivity behavior exactly,
  silence included. Corpora without near-baseline documents are unaffected.

* The word-level indices are unchanged: the `V+` restriction and the
  Section 3.8 support recommendations already guard that path. The gain
  table needs no change either; its paired gains were already
  complete-case, and the man page now states that documents excluded by the
  floor at either topic count drop from that pair.

# OpTop 0.14.0

The inference release: the held-out layer (Section 3.7 of Lewis and Grossetti 2026)
and the moment-based specification tests (Section 4) join the package, completing the
theory of the revised paper. A formula-level audit confirmed that the in-sample
surface of 0.13.0 conforms to the paper's July revision exactly.

### New Features

* **Held-out indices with conditional inference** (`optop_index_holdout()`): topic
  weights are folded in for independent evaluation documents, the harmonized support
  is rebuilt on the evaluation corpus over the augmented grid (grid plus the training
  baseline), and each index is reported with per-document scores, the Macro mean with
  its confidence interval, the Micro index and the Micro-Macro gap with delta-method
  standard errors, an optional stabilized score, and the paper's out-of-support
  convention (evaluation-only tokens enter the observed min-bin; document lengths
  keep every token). Zero-token evaluation documents are dropped with a warning.

* **Epsilon-adequacy selection** (`optop_gain_table()`): paired adjacent held-out
  gains on the actual grid with one-sided upper confidence bounds and the selection
  of the smallest K whose gain bound falls below the tolerance, with the documented
  post-selection caveats.

* **Moment-based specification tests** (`optop_moment_test()`): held-out residual
  probability vectors projected onto training-built, exactly mean-zero vocabulary
  instruments; three families (frequency contrast, frequency strata against the
  highest-frequency reference, fit-stratified via the training word-level deviance
  index) with Wald chi-square statistics, the scalar t equivalent, and marginal
  per-stratum t tests with optional Bonferroni or Benjamini-Hochberg adjustment.
  The moments are computed without forming any documents-by-vocabulary temporary.

* **text2vec WarpLDA support**: `optop_warplda()` wraps a WarpLDA fit with the
  document-topic matrix the sampler returns from `fit_transform()`, making the
  paper's estimator a first-class engine across the package; raw R6 objects fail
  with a pointer to the wrapper instead of falling through to the topicmodels
  adapter. The obsolete `LDA_t2v` stub is gone.

* **Fold-in adapters** (`optop_fold_in()`, internal): topicmodels via
  `posterior(newdata)`, seededlda via refitting with the trained model, WarpLDA via
  `$transform()`, NLPstudio via the stored backend fit.

* `optop_make_baseline()` gains additive smoothing (`smooth_lambda`) for the
  held-out support convention; `optop_make_partition()` accepts an external
  `pi_glob`, so held-out partitions use the training baseline as the null member.

### Fixes

* **Runtime and dependencies**: OpTop now requires R (>= 4.2.0), matching the
  package's current vignette tooling and shipped serialized data. The unused
  `tibble` import and obsolete `LazyData` declaration have been removed.

* The document kernel accumulates non-rare cells directly instead of as full minus
  rare sums: a rare cell with a large count on an eps-floored expectation (the
  held-out out-of-support column is exactly that case) made the subtraction cancel
  catastrophically in the Pearson family. In-sample results change at most in the
  last floating-point digits.

### Documentation

* Fitted probabilities follow the paper's notation (`p^K_jw`); the partition page
  states the augmented-grid union explicitly. The audit trail and the shrunken
  roadmap (cross-fitting helpers) are recorded in the development audit file.

* The man pages of the held-out family (`optop_index_holdout()`,
  `optop_gain_table()`, `optop_moment_test()`, `optop_warplda()`, and the internal
  `optop_fold_in()`) are written for a broad audience: each opens with the question
  the tool answers, walks through the pipeline step by step, maps every returned
  element to one sentence, and ships a complete runnable example.

* The vignette gains the section "Held-Out Validation and Topic-Number Choice": a
  seeded training and evaluation split of the inaugural corpus, a fresh
  training-split grid, the Macro held-out index with its confidence band, the
  epsilon-adequacy selection, and a frequency-strata moment test with a guided
  reading of the output.

# OpTop 0.13.0

The alignment release: the goodness-of-fit indices, the harmonized partition, and
the baseline now follow Lewis and Grossetti (2026), "Goodness-of-Fit Indices and
Diagnostics for Topic Models" (working paper), exactly.

### Breaking Changes

* **The harmonized rare-word set includes the no-topics baseline**: a word is now
  rare in document `j` when `min(pi_glob(w), min_K i_jw) < tau_j`, as the paper
  defines `C*_j`. All non-min expected counts are therefore bounded below by `c` on
  the active support, for the fitted models and the baseline alike. Masks, and hence
  all index values, change relative to 0.12.0.

* **Pearson min-bin inclusion rule**: the collapsed min-bin enters the chi-square
  discrepancy only when `min(min_K E_min, B_min) >= c`, decided once for the whole
  grid by `optop_make_partition()` (new `chisq_min_ok` field) and applied to the
  fitted and the null side simultaneously. The share of documents affected and the
  excluded observed mass are reported. Partitions computed by earlier versions keep
  the old always-include behavior with a warning.

* **Word-level deviance is the Poisson unit deviance**: the fitted side gains the
  linear correction `-(N - E)`, making every summand nonnegative and bounding the
  word-level index by one. The null side is unchanged, since the correction vanishes
  identically for the corpus baseline.

* **Aggregation restricts to non-degenerate units**: the Micro index sums over the
  documents (words) with positive baseline discrepancy only, and degenerate units
  now carry `NA` instead of `0` in `r2_doc` and `r2_word`.

* **Default `c` is now 1** in `optop_make_partition()` and `optop_index_table()`,
  the value the paper adopts and recommends when the deviance index is primary;
  `c = 5` remains appropriate for Pearson-primary work and explicit calls are
  unaffected.

### Deprecations (removal before v1.0.0)

* `ztest`: the methodology reserves inference for held-out evaluation; the
  in-sample Z-test is deprecated and full-corpus indices are documented as
  descriptive.
* `reopt` and `add_baseline_topic`: both enforce non-negativity, while the
  methodology treats negative indices as informative. `add_baseline_topic`
  defaults to `FALSE` (numerically a no-op under `reopt = "none"`).

### Minor Changes

* Document- and word-level results gain `d_model` and `d_null`, the per-unit fitted
  and baseline discrepancies, supporting the paper's recommended fit-versus-baseline
  scatter diagnostics.
* The index documentation is rewritten in the paper's terms: grid dependence of the
  harmonized support, per-family interpretation bounds, Micro as the
  discrepancy-weighted average, and the Micro-Macro gap as a covariance diagnostic.
* The efficiency contract of 0.12.0 is preserved: the baseline enters the partition
  kernel as the seed of the running minimum (no extra cost), the inclusion rule adds
  one blocked gemm sweep per model at partition time only, and all kernels remain
  bit-identical for any `n_threads`.

# OpTop 0.12.0

The C++ performance release: both hot paths of `optimal_topic()` now run in
compiled code under OpenMP.

### Major Changes

* **Parallel statistic core**: the per-document sort/cumsum/Pearson loop — the
  profiled bottleneck since the blocked-gemm rewrite of 0.9.7 — runs under OpenMP
  within each document block. Per-document results are reduced in a fixed order, so
  the output is bit-identical for any thread count.

* **The bootstrap moves to C++** (`calibrate = "bootstrap"`): the multinomial sampling
  is fused with the Pearson reduction (no `bins × B` count and deviance temporaries)
  and parallelized over documents. Each document owns a
  private RNG stream seeded deterministically from `(seed, document index)`, so results
  are bit-identical for any thread count and reproducible under `seed`; `set.seed()`
  at the R level still governs when `seed = NULL`.

* **New `n_threads` argument** on `optimal_topic()` (default `1L`), forwarded to both
  cores. It only affects wall time, never the numbers. On toolchains without OpenMP
  (e.g. the default macOS compiler) the cores run single-threaded; the thread-free
  C++ bootstrap speedup is kept.

* **Adaptive bootstrap sampler**: each document picks its sampling strategy from the
  data alone. Wide envelopes (more bins than tokens, the common regime on large
  vocabularies) draw the tokens directly through a Walker alias table at O(N_j) per
  replicate, with a closed-form Pearson contribution for untouched bins; narrow
  envelopes keep the conditional-binomial method with an early exit once the count is
  exhausted. Thread-count invariance and seed reproducibility are unaffected.

* **Envelope selection replaces the full sort**: the statistic locates the top-`P_j`
  head by `nth_element` selection with incremental slice sorts instead of sorting all
  `W` fitted probabilities per document.

* **Deterministic tie resolution**: envelope order is pinned to (probability
  descending, index ascending), the semantics of R's stable `order()`. Results are now
  identical across platforms even when fitted probabilities contain exact ties (as
  Gibbs estimates routinely do); previously the boundary of the envelope could differ
  between toolchains.

* **RNG stream change**: the compiled bootstrap draws through its own generator
  (splitmix64-seeded `std::mt19937_64`) instead of R's `rmultinom()`.
  Bootstrap-calibrated p-values therefore agree with 0.11.0 up to Monte-Carlo noise
  (order `1/sqrt(n_boot)`), not bit for bit. The test suite checks the new stream
  against the pure-R oracle in distribution and against the exact Haldane moments.

* **The discrepancy indices move to C++**: `optop_index_se()`, `optop_index_chisq()`,
  `optop_index_deviance()` and `optop_index_table()` are evaluated by a fused OpenMP
  engine (`src/index_core.cpp`) that computes every requested metric, for both the
  model and the no-topics baseline side, in a single traversal of each document block,
  at the document and at the word level. `optop_index_table()` now performs one sweep
  per model for the whole metric set, and the baseline side is still computed once per
  grid. The support decomposition and the epsilon conventions are reproduced exactly;
  results are unchanged and bit-identical for any thread count.

* **Compiled harmonized partition**: `optop_make_partition()` computes the running
  minimum across the model grid one vocabulary block at a time in compiled code
  (`src/partition_core.cpp`). The `J × W` minimum matrix is no longer materialized:
  peak memory is the returned rare mask plus one `J × block` buffer.

* **`n_threads` on the index family**: `optop_index_se()`, `optop_index_chisq()`,
  `optop_index_deviance()`, `optop_index_table()` and `optop_make_partition()` gain an
  `n_threads` argument (default `1L`) with the same contract as `optimal_topic()`:
  wall time only, never the numbers.

* **Bootstrap envelope marshalling**: `optimal_topic()` passes the statistic core's
  flattened envelope export directly to the bootstrap core; the per-document list
  round trip and its transient copy (about 1 GB on the largest benchmark corpus) are
  gone. The list split remains only in the moment-matching path, which needs it.

### Minor Changes

* New vignette section **"A Note On Computational Efficiency"** reporting a simulation
  study (speed, thread scaling up to the machine's core count, peak memory) across four
  synthetic corpus scales, up to J = 10,000 documents and W = 20,000 features, with
  the pre-0.12.0 R bootstrap as baseline; reproducible via
  `data-raw/benchmark-efficiency.R`. The study also covers the harmonized partition
  and the three-metric index table, with the pre-engine R implementation as baseline.

* `src/Makevars` and `src/Makevars.win` add `$(SHLIB_OPENMP_CXXFLAGS)`.

* `optimal_topic()` keeps the adapter extractions for the evaluation loop under a
  fixed memory budget instead of extracting every model twice.

* New tests: bit-identical thread sweeps for all compiled cores (statistic, bootstrap,
  index engine, partition kernel), seed reproducibility of the compiled bootstrap,
  distributional agreement with the `rmultinom()` oracle, and an exact tie-resolution
  test against the R `order()` reference.

# OpTop 0.11.0

### New Features — engine generalization

* **`optimal_topic()` accepts any supported topic model, not just `LDA_VEM`**: the first
  argument is now `topic_models` and takes, through the `optop_as_theta_phi()` adapter
  family, `topicmodels::LDA()` fits (VEM **and Gibbs**), `topicmodels::CTM()` fits, all
  three seededlda models (`textmodel_lda()`, `textmodel_seededlda()`,
  `textmodel_seqlda()`) and NLPstudio fits (`nlp_topic_fit`, via the stored `dtw`/`tww`
  weights or, when absent, the raw backend `model_object`). Engines can be mixed within
  one grid fitted on the same corpus. No statistical change is involved: the Test 1
  statistic consumes each model only through its fitted word probabilities Θ Φ.

* **The C++ core is model-agnostic** (and simpler): it receives plain `theta`/`phi`
  matrices — no S4 slot access in compiled code — and the weighted dfm transposed once
  for the whole grid. The hot path is now free of R API calls, the precondition for the
  planned OpenMP parallelization. The pre-0.9.9 pipeline behind `selection = "legacy"`
  is frozen verbatim in its own translation unit (`optimal_topic_core_legacy.cpp`),
  remains VEM-only, and will be deleted together with the rule before v1.0.0.

### Stronger input contracts

* **Vocabulary alignment is now validated**: all models of a grid must share one
  vocabulary in one order; `weighted_dfm` features that are the same set in a different
  order are reordered automatically (signalled), any other mismatch is an error.
  Previously a feature-permuted dfm was silently mis-scored.

* **Document alignment is identifier-based, per model**: every model must expose unique
  document identifiers (a mandatory field of the adapter contract — positional
  alignment is never assumed). Documents missing from *any* model are dropped with a
  warning (previously only the first model was checked), and each retained dfm row is
  matched to that model's own theta row, so engines that saw the documents in different
  orders are scored correctly.

* Grids are sorted by increasing K (signalled) and duplicated K values are an error; a
  `weighted_dfm` whose rows do not sum to 1 is flagged as well.

### Deprecations

* `optimal_topic(lda_models = )` is deprecated in favor of `topic_models`; the old name
  keeps working with a lifecycle warning and will be removed before v1.0.0.

### Internals

* The adapter contract is extended to `list(theta, phi, K, docs, terms)` and enforced
  by a dedicated validator (simplex rows, unique identifiers, consistent dimensions);
  unsupported classes now fail with the list of supported engines.

* seededlda joins Suggests (adapter tests cover all three constructors); the NLPstudio
  adapter is tested against a hand-built object, so OpTop takes no dependency on
  NLPstudio.

# OpTop 0.10.1

### Minor Changes

* **Internal helpers fully documented**: the calibration engine, the index workers, the
  alignment validators and the `optop_as_theta_phi()` adapter family now carry proper
  roxygen documentation with `@keywords internal` (browsable via `?OpTop:::name`),
  including the calibration mathematics typeset with `\eqn{}`.

* Remaining documentation mathematics (the cross-document Z-test) moved from unicode
  prose to `\eqn{}`; references extended — Satterthwaite (1946) and Davison & Hinkley
  (1997) join `?optimal_topic`, and the index/partition pages now cite Lewis & Grossetti
  (2022) and Agresti (1996).

* Dead code removed: ten stale `utils::globalVariables()` declarations referenced by
  no live or deprecated code, plus a dangling roxygen fragment in `R/utils.R`.

* **Package overview rewritten**: `?OpTop-package` now presents the current toolkit
  (Test 1 with its selection rules and calibration, the discrepancy indices, the
  harmonized partition) and cites Lewis and Grossetti (2022) — the previous page still
  referenced a 2019 preprint and predated the index family.

# OpTop 0.10.0

### New Features — calibrated p-values for Test 1

* **`optimal_topic()` gains a `calibrate` argument** replacing the chi-square yardstick
  with the distribution of the statistic under the *fitted-model null* (each document
  multinomial of its own length with the model's fitted word probabilities, Θ and Φ held
  fixed):
  - `calibrate = "bootstrap"`: parametric bootstrap drawn exactly on the collapsed
    envelope bins (`n_boot` replicates, optional `seed`), empirical upper-tail p-value
    with the `(1 + #{T* >= T}) / (n_boot + 1)` correction. With calibration, `alpha` is a
    genuine Type-I error rate under the conditional null.
  - `calibrate = "moment"`: exact Haldane (1937) multinomial moments matched to a
    Satterthwaite scaled chi-square — closed form, no simulation.
  - New `doc_lengths` argument (e.g. `quanteda::ntoken()` of the counts dfm) supplies the
    document lengths the null depends on; `pval` carries the calibrated value the
    selection rules use, with the asymptotic `pval_chisq` kept alongside.
  - The C++ core exports the per-document envelope (additive `return_envelope` flag; the
    hot path is untouched), and the bootstrap costs about `n_boot x df` vectorized flops
    per model. The man page documents the full statistical reasoning, including the
    fixed-Θ̂Φ̂ caveat.

### Major Changes

* **`verbose = TRUE` is the new default** of `optimal_topic()` (reversing the 0.9.8
  silent default): the cli report — setup, calibration announcement, progress bar,
  selection summary, timing — is now shown unless `verbose = FALSE`.

# OpTop 0.9.10

### New Features

* **Vignette**: `vignette("OpTop")` walks through optimal-K selection on the U.S.
  Presidential Inaugural Address corpus (both `selection = "sequential"` and
  `selection = "min"`, and the role of `q`), the discrepancy indices, and a simulation
  with known true K. The heavy LDA grid is precomputed by `data-raw/vignette-data.R`;
  only the small derived tables ship with the package.

### Minor Changes

* `optimal_topic()`'s plot now uses `ggplot2::theme_bw()` and carries a subtitle with
  the selection method and the selected K.

* Style: no spaces inside parentheses in the maintained files.

* `AUDIT.md` records the proposal for calibrating Test 1's null distribution
  (parametric-bootstrap null, moment matching, count-based statistic), planned as a
  dedicated future release.

# OpTop 0.9.9

### Major Changes — `optimal_topic()` calibrated to the published test

**This release changes numerical results.** The implementation is now faithful to
Test 1 (Equation 8) of Lewis and Grossetti (2022):

* **Statistic**: each document's Pearson term is scaled by `P_j + 1` (the number of bins,
  previously `P_j`), and the corpus statistic is chi-square with `df = sum_j P_j`. The
  returned `OpTop` column reports the standardized statistic (raw / `df`, the scale of the
  paper's Figure 2) and the table gains a `df` column.

* **P-values**: upper tail of the raw statistic on its full degrees of freedom
  (previously: lower tail of the standardized statistic on 1 df, which asked "is the fit
  suspiciously perfect?" rather than "does the model fit?").

* **Selection rules**: new `selection` argument — `"sequential"` (default; smallest `K` the
  test fails to reject at `alpha`, with a warned fallback to the global minimum),
  `"min"` (the global minimum of the standardized statistic — the rule used in the paper's
  case study), and `"legacy"` (the frozen pre-0.9.9 pipeline, bit-identical to v0.9.8).
  `"legacy"` is deprecated at birth and will be removed before v1.0.0.

* **Cutoff**: the "relatively important" words are now the smallest sorted head whose
  cumulative mass strictly exceeds `q`, keeping the crossing word, so the collapsed tail
  stays strictly below `1 - q` (the paper's footnote 5); the 4-decimal rounding hack is
  gone. The default `q` moves from `0.80` to `0.95`, matching the paper's numerical
  setting `I^K = 0.05`.

* Documentation now includes an honest calibration caveat: with `df = sum_j P_j` in the
  thousands the chi-square p-values saturate near 0/1 unless the fit is borderline, in
  which case the shape and minimum of the standardized curve carry the information.

# OpTop 0.9.8

### Major Changes

* **`optimal_topic()` is silent by default and reports via `cli`**: the new `verbose`
  argument (default `FALSE`) replaces the old always-on console chatter. With
  `verbose = TRUE` the function shows a `cli` header, a setup summary, a live progress
  bar across the model grid (one bar instead of one line per model — relevant for large
  grids such as K = 10..500), the selected model with the selection rule, and the elapsed
  time. Dropping documents that the models never saw is always signalled, regardless of
  `verbose`. The C++ core no longer prints anything.

### Deprecations

* `topic_stability()`, `agg_topic_stability()`, `agg_document_stability()` and
  `get_topic_models()` are deprecated (with `lifecycle` warnings and badges) and
  scheduled for removal in a future release, together with the internal `topic_match()`.
  The package is converging on the discrepancy-index API (`optop_index_*()`).

### Minor Changes

* Style modernization in the optimal-topic and discrepancy files: `<-` for assignment
  and fully qualified `pkg::fun()` calls instead of `@importFrom` (the `data.table`
  import stays for its non-standard evaluation).

* `cli` joins the imports; `lifecycle` moves from Suggests to Imports.

* README: release badge updated, `verbose` shown in the Quick Start, and the "Current
  support and extensions" section rewritten to reflect the actual implementation
  (`optimal_topic()` requires VEM fits; the discrepancy indices accept VEM and Gibbs
  fits via `optop_as_theta_phi()` adapters; further adapters are planned).

* `DESCRIPTION` now cites the methodological reference, Lewis and Grossetti (2022)
  <https://jmlr.org/papers/v23/19-297.html>, following CRAN guidelines.

* Continuous integration: GitHub Actions workflows for `R CMD check` (Windows, macOS
  and Ubuntu across oldrel/release/devel) and test coverage via `covr`/Codecov, with
  the corresponding badges in the README.

# OpTop 0.9.7

### Major Changes

* **`optimal_topic()` core rewritten around blocked matrix products**: the per-document
  hot loop used to build a dense $W \times K$ temporary and row-sum it; that row sum is one
  row of $\gamma \, e^{\beta}$, so each block of documents now needs a single BLAS product
  per model. Documents are densified from a transposed copy of the dfm in contiguous
  column blocks instead of element-by-element sparse row reads (which were also redone
  for every model). Results are unchanged; see `AUDIT.md` for benchmarks.

### Bug Fixes

* **Document alignment in `optimal_topic()`**: only *membership* of the dfm documents in
  the models' `@documents` was checked, never *order*, so a dfm whose rows were permuted
  relative to the fitted models was silently mis-scored. Each dfm row is now paired with
  its `@gamma` row by document identifier.

* **Plotting warnings removed**: `optimal_topic()` and `topic_stability()` plots no longer
  trigger the ggplot2 `size`-for-lines deprecation (now `linewidth`) nor the length-1
  aesthetics warning (the optimum marker is drawn with `annotate()`). The ggplot2
  requirement is now `>= 3.4.0`.

### Minor Changes

* Extended the `testthat` suite to `optimal_topic()`, with a naive per-document reference
  implementation mirroring the C++ semantics, document-permutation and document-removal
  invariance checks, and input validation tests.

* Re-documented with `roxygen2` 8.0.0 and fixed a typo in the maintainer's email address
  in `DESCRIPTION`.

* Added `AUDIT.md` tracking the audit of the optimal-topic pipeline: what is fixed, what
  is methodological and deliberately unchanged (lower-tail p-value, rounded cumulative-mass
  cutoff), and what is deferred (OpenMP, model-class generalizability).

# OpTop 0.9.6

### New Functions

* **Discrepancy-based goodness-of-fit indices**: Four new functions implement regression-style 
  $R^2$ indices for topic model evaluation:
  - `optop_index_se()`: Squared-error index.
  - `optop_index_chisq()`: Pearson chi-square index.
  - `optop_index_deviance()`: Deviance index.
  - `optop_index_table()`: Computes indices across a grid of models and returns a data.table.

### Major Changes

* **Word-level aggregation**: All index functions support `level = "word"` to compute per-word 
  R² indices with Micro-Word (frequency-weighted) and Macro-Word (unweighted) aggregations. 
  This reveals which words are well-captured vs. poorly modeled by the topic structure.

* **Z-test for cross-document inference**: Setting `ztest = TRUE` tests whether the topic model 
  provides statistically significant improvement over the no-topics baseline.

* **Block-based processing**: All index functions now use memory-efficient block-based computation 
  for both word-level and document-level aggregation. This enables processing of large corpora 
  (100K+ documents) without memory issues. The `block_size` parameter controls word-level blocking;
  document-level blocking is adaptive and computed internally. The `block_size` parameter is adaptive.

* **Vectorized document-level computation**: Replaced per-document for loops with vectorized 
  matrix operations using the decomposition $D = D_{full} - D_{rare} + D_{min}$, significantly 
  improving performance on large corpora.


* **Baseline discrepancy computed once per grid**: `optop_index_table()` now computes the 
  model-independent no-topics discrepancy $D_{null}$ once per metric and reuses it across the 
  whole grid of $K$, roughly halving the per-model cost on large corpora. Results are unchanged.

### Minor Changes

* Added a `testthat` suite covering the partition, baseline, and all index functions (document 
  and word level, macro, Z-test, re-optimization, block-size invariance, alignment guards), 
  including an independent naive reference implementation of the indices that follows the 
  paper's definitions with explicit loops.

* Vectorized the word-level deviance computation and removed redundant temporaries in the 
  document-level blocks; numerically identical results.

* `optop_index_table()` now returns its `data.table` visibly, so the result prints at the console.

* Updated documentation to Markdown format with better rendering of mathematical equations.

* Fixed CRAN notes and warnings.

* Removed unused internal helpers and dead code in `optop_make_partition()`.

### Bug Fixes

* Fixed S4 subscript error when using `quanteda::dfm` or `Matrix::dgCMatrix` inputs. The issue 
  occurred because extracted values retained S4 class; now handled with explicit `as.numeric()` 
  conversion.

---


# OpTop 0.9.5

### Major Changes

* `word_proportions()`: after careful testing, we decided to deprecate this function. 
The main reason is the lack of support for high-dimensional `dfm` inputs. Now all the
functions that require word proportions take advantage of the class `dfm` from **quanteda**.

* `sim_dfm()`: we introduced a new function to easily simulate a document-feature-matrix from
a LDA specification. This is useful for simulating corpora for testing.

### Minor Changes

* Document matching between input LDAs and `weighted_dfm` has been improved. Now the check does not
rely on a specific `docvar` anymore but rather on the internal document naming convention as
defined in the `corpus`.

* Small tweaks in the computation speed for certain operations.

* Improved documentation.

### Bug Fix

* If the optimal model is the last in the list of LDA objects passed to `optimal_topic()`, now all the 
other functions devoted to testing return `NULL` with a message. This is because there is nothing
to test in terms of topic stability above the optimal model.


---

# OpTop 0.9.4

### Major Changes

* `word_proportions()` is much general now and can carry out complex 
preprocessing routines. 
* `optimal_topic()` is much faster now due to matrix implementation. Now it checks
for the presence of documents in both the `corpus`/`dfm` and the estimates by `LDA()`.

### Minor Changes

* The argument `remove_documents` in `word_proportions()` is now set to `FALSE`
by default. This automatically triggers document check in `optimal_topic()`. 
* Preprocessing in `word_proportions()` is achieved by the new general argument
`...`.
* `word_proportion()` can now work on a `corpus` and `dfm` objects as defined
in __quanteda__.
* Improved documentation 

---

# OpTop 0.9.3

* Fixing tiny minor things

---

# OpTop 0.9.1

* Just kidding! No more data here...we'll figure out something later on...

---

# OpTop 0.9.0

### Major Changes

* A set of pre-run LDA models are now available.

# OpTop 0.8.0

### Major Changes

* `agg_document_stability()` has been fully implemented. The function returns
both the Aggregated Document Stability test and the F-test on informative
and uninformative components.

### Minor Changes

* General improvements in plots.

* Better documentation.

* Slighty faster functions.

---

# OpTop 0.7.0

### Minor Changes

* `agg_topic_stability()` can now compute smoothed tests and plot the results
accordingly.

* Support for final convertion to a `tibble` table spreaded out to all
functions,

* All eligible functions get better plots.


---

# OpTop 0.6.0

### Major Changes

* `optimal_topic()` gains the parameter `q` which allows to select the quantile
of the cumulative probability of word weights to consider as relevant.

* `optimal_topic()` now finds the optimal number of topics either by significance
levels or by forcing the algorithm to reach the global minimum. This is 
controlled by the new parameter `alpha`. 

### Minor Changes

* `optimal_topic()` drops both `threshold` and `q_type`.

* In `optimal_topic()`, `convert` now supports `tibble` structure.

---

# OpTop 0.5.0

### Major Changes

* Function `agg_topic_stability()` has been widely improved.

* All functions which return a test now gain the new argument `do_plot`. This
plot the test statistic as a function of the number of topics.

* The argument `test` has been removed from `topic_stability()` which now 
returns only the aggregate statistic. 

* The argument `compute_res` has been finally removed from `topic_stability()`.

* `topic_stability()` now returns either a data.frame or a data.table with 
the LDA specifiction associated to each statistic (i.e. column `topic`).

### Other Changes

* Improved documentation for some functions.

---

# OpTop 0.4.0

### New Functions

Since we have two more functions, I feel like this deserves a jump in 
package version. 

* `topic_match()`: detect and extract informative and uninformative components.

* `agg_topic_stability()`: implements _Test 4_ from the methodological paper 
[Lewis and Grossetti (2019)].

### Other Changes

* Improved documentation for some functions.

---

# OpTop 0.2.0

### New Functions

* `get_topic_models()`: handy function to immediately get the list of topic models
the user wants to process from a specified environment;

* `topic_stability()`: implements _Tests 2_ and _3_ from the methodological paper 
[Lewis and Grossetti (2019)].

### Other Changes

* Formal declaration of `LDA_VEM` objects as functions input.

* All the functions now have more detailed and better documentations.

* Added Continuous Integration with Travis CI and AppVeyor.

### Bug Fixes

* Choice of quantile algorithms is now fully supported.

---

# OpTop 0.1.0

First version! 

### New Functions

* `word_proportions()`: computes word proportions from a `corpus` object created 
by __quanteda__ [Benoit et al. (2018)];

* `optimal_topic()`: implements _Test 1_ from the methodological paper 
[Lewis and Grossetti (2019)].
