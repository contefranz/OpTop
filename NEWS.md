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
