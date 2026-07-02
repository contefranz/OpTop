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

### Confirmed, deliberately unchanged (methodological — author's call)

- **`pchisq(stat, df = 1, lower.tail = TRUE)`.** A goodness-of-fit p-value is
  conventionally the upper tail. With the lower tail, a *small* statistic
  yields a *small* p-value, which interacts with the `alpha` selection rule
  ("first model with p-value ≤ alpha"). This is tied to the paper's selection
  procedure and is left exactly as is; the test suite encodes the current
  behavior.
- **`round(val * 1e4) / 1e4 > q` cumulative-mass cutoff.** Rounding to 4
  decimals before comparing with `q` makes the truncation index unstable near
  ties (documents whose cumulative mass sits within 5e-5 of `q` can flip).
  Kept verbatim to preserve numerical output.

## Deferred

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

## Benchmark

_To be filled in: before/after timings of `optimal_topic()` on a simulated
corpus (J ≈ 2000, W ≈ 5000, grid of K), same machine, `do_plot = FALSE`._
