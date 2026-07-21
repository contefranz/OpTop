#' OpTop goodness-of-fit indices for a single topic model
#'
#' Compute the goodness-of-fit indices \eqn{R^2_D(K)} of Lewis and Grossetti
#' (2026) for one fitted topic model: the proportional reduction in
#' discrepancy achieved by the \eqn{K}-topic model relative to the no-topics
#' baseline in which every document follows the corpus word distribution,
#' \deqn{R^2_D(K) = 1 - D(K) / D(\mathrm{null}),}{R2_D(K) = 1 - D(K)/D(null),}
#' for the squared-error, Pearson chi-square, and deviance discrepancy
#' families. Document-level indices are evaluated on the harmonized support
#' of [`optop_make_partition()`]; word-level indices are evaluated on the
#' unbinned vocabulary.
#'
#' @param model A single fitted `topicmodels::LDA` model (VEM or Gibbs), or
#'   a loader function returning one (materialized on demand).
#' @param dtm A counts document–term matrix (documents × terms). Use `Matrix::dgCMatrix`
#'   or a `quanteda::dfm` that represents raw counts (not weighted), or an
#'   [optop_corpus()] of count shards for corpora past the size of a single
#'   matrix (shards stream one at a time; per-document results are
#'   bit-identical to an unsharded run).
#' @param partition Result from [`optop_make_partition()`], computed on a DTM that is
#'   *aligned* to the model vocabulary (same features and order).
#' @param baseline Result from [`optop_make_baseline()`], computed on the same aligned
#'   counts DTM used for `partition`.
#' @param macro Logical; if `TRUE`, also compute the Macro index (equal-weight average of
#'   the per-document indices). Default: `FALSE`.
#' @param reopt `r lifecycle::badge("deprecated")` Character scalar indicating
#'   optional re-optimization. The blending is not part of the methodology
#'   and is scheduled for removal before v1.0.0; use the default `"none"`.
#' @param level Character; aggregation level for the index. `"document"` (default) computes
#'   document-level indices aggregated across words. `"word"` computes word-level indices
#'   aggregated across documents.
#' @param block_size Integer or `NULL`; number of vocabulary terms to process at once
#'   during word-level computation. Smaller values use less memory but may be slower.
#'   If `NULL` (default), automatically chosen based on corpus size to target ~1.5 GB
#'   memory usage. Only used when `level = "word"`.
#' @param ztest `r lifecycle::badge("deprecated")` Logical; if `TRUE`, append
#'   the in-sample Z-test. The methodology reserves inference for held-out
#'   evaluation, so the in-sample test is scheduled for removal before
#'   v1.0.0; full-corpus indices are descriptive.
#' @param n_threads Integer; number of OpenMP threads used by the compiled
#'   index engine (default `1L`). Results are identical for any value; only
#'   wall time changes.
#' @param min_null Nonnegative scalar; the null-discrepancy floor of the
#'   document-level aggregation. Documents with baseline discrepancy
#'   \eqn{D_j(\mathrm{null}) < \mathrm{min\_null}}{D_j(null) < min_null}
#'   are excluded from the Micro and Macro aggregations and carry an `NA`
#'   document-level index. The default `NULL` resolves to the partition
#'   threshold constant `c`; `min_null = 0` restores the pre-0.14.1 strict
#'   positivity rule. See the *Null-discrepancy floor* paragraph in Details.
#'
#' @details
#' **Harmonized support.** For each document, rare words (as determined by
#' [`optop_make_partition()`], whose union includes the no-topics baseline)
#' are collapsed into a single "min" bin. Document-level indices are then
#' evaluated on `{non-rare terms} ∪ {min}`, so that changes in the index
#' across \eqn{K} reflect changes in model fit, not changes in binning.
#' Because the harmonized set is a union over the whole model grid, the
#' reported value at any fixed \eqn{K} depends on the pair (grid, `c`):
#' compare indices across studies only under a common grid and `c`.
#'
#' **Interpretation.** The indices are descriptive effect sizes anchored by
#' the null baseline (lower reference) and the saturated model (upper
#' reference), not tests of correct specification. The deviance index is
#' bounded above by one; the Pearson and squared-error indices are not
#' confined to \eqn{[0, 1]}, and negative values are informative: the fitted
#' model performs worse than the no-topics baseline under that discrepancy.
#' In-sample values computed on the estimation corpus summarize descriptive
#' fit and should not be the sole basis for choosing \eqn{K}.
#'
#' **Pearson min-bin rule.** For the chi-square family the collapsed min-bin
#' of a document enters only when
#' \eqn{\min(\min_K E^K_{j,\min}, B_{j,\min}) \ge c}{min(min_K E_min, B_min) >= c},
#' decided once for the whole grid by [`optop_make_partition()`] and applied
#' to the fitted and the null discrepancy simultaneously; excluded shares
#' are reported.
#'
#' **Alignment requirements.** The following must share the same vocabulary and column
#' order: the `dtm` passed to the index, the DTM used for `partition` and `baseline`, and the
#' model's term–topic matrix. If they differ, align with `optop_align_dtm_to_models()` and
#' recompute `partition` and `baseline`.
#'
#' **Counts only.** SE/chi-square/deviance indices are defined for multinomial counts.
#' Do not pass weighted matrices (e.g., proportions from `quanteda::dfm_weight(scheme = "prop")`
#' or tf-idf). If your workflow uses proportions elsewhere, reconstruct counts before calling.
#'
#' **Document-level aggregation** (`level = "document"`). The Micro index
#' pools discrepancies before the ratio and is therefore a
#' discrepancy-weighted average of the document-level fits, with weights
#' proportional to the baseline discrepancy \eqn{D_j(\mathrm{null})}{D_j(null)};
#' the Macro index is their unweighted average. Both aggregate over the
#' non-degenerate documents
#' \eqn{J_+ = \{j : D_j(\mathrm{null}) \ge \mathrm{min\_null}\}}{J+ = {j: D_j(null) >= min_null}}
#' only, and excluded documents carry an `NA` document-level index. The
#' Micro–Macro gap is itself diagnostic: it equals the covariance between
#' document-level fit and baseline discrepancy (normalized by the mean
#' baseline discrepancy), so a positive gap indicates that fit concentrates
#' in high-discrepancy (often long or atypical) documents. The returned
#' `d_model` and `d_null` vectors support the recommended scatter of
#' document-level fit against baseline discrepancy or length.
#'
#' **Null-discrepancy floor.** At partition resolution
#' \eqn{\tau_j = c/L_j}{tau_j = c/L_j} a document whose support is entirely
#' swept into the collapsed min-bin has observed equal to baseline counts by
#' construction, so its baseline discrepancy is zero in every family at
#' once, and a document whose word distribution sits close to the corpus
#' baseline has \eqn{D_j(\mathrm{null})}{D_j(null)} arbitrarily near zero.
#' In both regimes the document-level index
#' \eqn{1 - D_j(K)/D_j(\mathrm{null})}{1 - D_j(K)/D_j(null)} is unbounded,
#' and the Macro index, an unweighted mean, can be destroyed by a handful of
#' such documents, while the Micro index is far less exposed because it
#' weights documents by \eqn{D_j(\mathrm{null})}{D_j(null)} itself. The
#' aggregation therefore conditions on
#' \eqn{D_j(\mathrm{null}) \ge \mathrm{min\_null}}{D_j(null) >= min_null},
#' with default equal to the partition constant `c`: the estimand becomes
#' the average fit among documents whose baseline discrepancy is at least
#' the resolution of the evaluation support. Tiny baseline discrepancies
#' signal collapsed supports rather than sampling noise: under the null, a
#' document's binned baseline discrepancy is approximately chi-square with
#' (number of effective bins - 1) degrees of freedom, so values orders of
#' magnitude below one arise only when the partition has degenerated.
#' Exclusions are reported once per call and returned as
#' `n_null_excluded`; report sensitivity to `min_null` alongside `c` when
#' the excluded share is non-negligible.
#'
#' **Word-level aggregation** (`level = "word"`). Word-level indices are
#' computed on the unbinned vocabulary and are descriptive diagnostics: they
#' reveal which words the topic structure captures well and which it
#' mispredicts. The word-level deviance is the Poisson unit deviance
#' \eqn{2\sum_j [N_{jw}\log(N_{jw}/E_{jw}) - (N_{jw}-E_{jw})]}{2 sum_j [N log(N/E) - (N - E)]},
#' which is nonnegative term by term, so the word deviance index never
#' exceeds one; for the corpus baseline the linear correction vanishes
#' identically. Micro-Word weights words by their baseline discrepancy;
#' Macro-Word averages the word-level indices over the words with positive
#' baseline discrepancy (degenerate words carry `NA`). Restrict attention to
#' words with adequate support when interpreting single-word values.
#'
#' @return When `level = "document"`, a list with:
#' - `r2`: scalar Micro index over \eqn{J_+}{J+}.
#' - `r2_macro`: scalar Macro index if `macro = TRUE`, otherwise `NULL`.
#' - `r2_doc`: numeric vector (length J) with per-document indices; `NA` for
#'   documents with zero baseline discrepancy.
#' - `d_model`, `d_null`: numeric vectors (length J) with the per-document
#'   fitted and baseline discrepancies.
#' - `n_null_excluded`, `null_excluded_share`: number and share of documents
#'   excluded by the null-discrepancy floor.
#' - `K`: number of topics in `model`.
#' - `metric`: one of `"se"`, `"chisq"`, `"deviance"`.
#' - `ztest`: (if `ztest = TRUE`; deprecated) list with `z`, `pval`, `se`, `ci`, `J`.
#'
#' When `level = "word"`, a list with:
#' - `r2_word`: named numeric vector (length W) of per-word \eqn{R^2_{D,w}(K)};
#'   `NA` for words with zero baseline discrepancy.
#' - `r2_micro_word`: scalar Micro-Word index (discrepancy-weighted average).
#' - `r2_macro_word`: scalar Macro-Word index (unweighted average over words
#'   with positive baseline discrepancy).
#' - `d_model`, `d_null`: numeric vectors (length W) with the per-word fitted
#'   and baseline discrepancies.
#' - `K`: number of topics in `model`.
#' - `metric`: one of `"se"`, `"chisq"`, `"deviance"`.
#'
#' @examples
#' \dontrun{
#' 
#' library(quanteda)
#' library(OpTop)
#' 
#' K_grid = 2:10
#' 
#' # Tokenize the corpus
#' toks <- data_corpus_inaugural %>% 
#'   tokens(remove_punct = TRUE, 
#'          remove_symbols = TRUE, 
#'          remove_numbers = TRUE) %>% 
#'   tokens_tolower() %>% 
#'   tokens_remove(stopwords())
#'   
#' # Create the document-feature-matrix
#' mydfm <- dfm(toks)
#' 
#' # Estimate topic models via VEM
#' VEM_models <- lapply(
#'   K_grid, function(k) {
#'     topicmodels::LDA(x = mydfm, k = k)
#'   }
#' )
#'
#' # Same DTM used at fit and evaluation --> no separate alignment here
#' partition <- optop_make_partition(models = VEM_models, dtm = mydfm, c = 1)
#' baseline  <- optop_make_baseline(dtm = mydfm)
#'
#' # Choose model with k = 10
#' # SE
#' res_se <- optop_index_se(
#'   model     = VEM_models[[9]],
#'   dtm       = mydfm,
#'   partition = partition,
#'   baseline  = baseline
#' )
#' res_se$r2
#'
#' # Chi-square
#' res_x2 <- optop_index_chisq(VEM_models[[9]], mydfm, partition, baseline)
#'
#' # Deviance
#' res_dev <- optop_index_deviance(VEM_models[[9]], mydfm, partition, baseline)
#' }
#'
#' @seealso
#' [`optop_make_partition()`], [`optop_make_baseline()`], [`optop_index_table()`];
#' internal helper: `optop_align_dtm_to_models()`.
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2026). Goodness-of-fit indices and
#' diagnostics for topic models. Working paper.
#'
#' Lewis, C. M. and Grossetti, F. (2022). A statistical approach for optimal
#' topic model identification. *Journal of Machine Learning Research*,
#' 23(58), 1--20. <https://jmlr.org/papers/v23/19-297.html>
#'
#' Agresti, A. (1996). *An Introduction to Categorical Data Analysis*.
#' Wiley, New York.
#'
#' @name optop_index
#' @aliases optop_index_se optop_index_chisq optop_index_deviance
NULL
#' @describeIn optop_index Squared-error index \eqn{R^2_{SE}}{R2_SE}.
#' @param add_baseline_topic `r lifecycle::badge("deprecated")` Logical; if
#'   `TRUE`, augment topics with a baseline row to force non-negativity when
#'   combined with `reopt = "se"`. The methodology treats negative indices as
#'   informative, so the augmentation is scheduled for removal before v1.0.0.
#'   Default: `FALSE` (no effect under `reopt = "none"`).
#'   (Only for `optop_index_se()` and `optop_index_chisq()`).
#' @export
optop_index_se <- function(model, dtm, partition, baseline,
                           macro = FALSE, reopt = c("none", "se"),
                           add_baseline_topic = FALSE,
                           level = c("document", "word"),
                           block_size = NULL,
                           ztest = FALSE,
                           n_threads = 1L,
                           min_null = NULL) {

  level <- match.arg(level)
  reopt <- match.arg(reopt)
  min_null <- .optop_resolve_min_null(min_null, partition$c)
  .optop_deprecate_index_args("optop_index_se", reopt, add_baseline_topic,
                              ztest)
  res <- .optop_index_se_impl(model, dtm, partition, baseline, macro, reopt,
                              add_baseline_topic, level, block_size, ztest,
                              n_threads = n_threads, min_null = min_null)
  if (level == "document") {
    .optop_report_null_floor(res$n_null_excluded, res$null_excluded_share,
                             min_null, "se")
  }
  res
}

#' Worker behind optop_index_se()
#'
#' Carries the actual squared-error computation; the exported wrapper only
#' validates and forwards. Shares the signature of [optop_index_se()] plus
#' `null_disc`.
#'
#' @param null_disc Optional model-independent baseline discrepancy
#'   (per-document \eqn{D_{null}}{D_null} when `level = "document"`, per-word
#'   SST when `level = "word"`), as computed by `.optop_index_null()`, so
#'   that grid evaluations across \eqn{K} reuse it instead of recomputing it
#'   for every model. Ignored when `reopt != "none"` because the SE
#'   re-optimization needs the full baseline counts per document.
#' @inheritParams optop_index_se
#'
#' @return The same list the exported wrapper documents.
#'
#' @keywords internal
.optop_index_se_impl <- function(model, dtm, partition, baseline, macro, reopt,
                                 add_baseline_topic, level, block_size, ztest,
                                 null_disc = NULL, n_threads = 1L,
                                 min_null = 0) {

  if (.optop_is_corpus(dtm) && identical(reopt, "se")) {
    stop(paste("the deprecated SE re-optimization is a small-corpus path",
               "and does not accept an optop_corpus"))
  }
  tp <- optop_as_theta_phi(.optop_materialize_model(model))
  vocab_model <- colnames(tp$phi)
  pi_row <- .optop_validate_alignment(vocab_model, dtm, partition, baseline)

  if (reopt != "none") null_disc <- NULL

  theta <- tp$theta
  phi <- tp$phi
  J <- nrow(theta)

  # optional baseline topic augmentation (guarantees R2 >= 0 for SE)
  if (add_baseline_topic) {
    phi <- rbind(phi, pi_row)
    # keep θ as-is; any θ-reopt may allocate weight to baseline topic
    theta <- cbind(theta, rep(0, J))
  }

  if (level == "document" && reopt == "se") {
    return(.optop_index_se_reopt(theta, phi, dtm, partition, pi_row,
                                 tp$K, macro, ztest, min_null))
  }

  eng <- .optop_index_engine(theta, phi, dtm, partition, pi_row, "se",
                             level, block_size, n_threads,
                             do_null = is.null(null_disc))
  if (level == "word") {
    d_null <- if (is.null(null_disc)) eng$se$null else null_disc
    return(.optop_index_result_word(eng$se$model, d_null, vocab_model,
                                    tp$K, "se"))
  }
  D_null <- if (is.null(null_disc)) eng$se$null else null_disc
  .optop_index_result_doc(eng$se$model, D_null, tp$K, "se", macro, ztest,
                          min_null)
}

#' SE re-optimization path (document level)
#'
#' The pre-engine block loop, kept verbatim: the per-document lambda blending
#' needs the full baseline counts, so it does not fit the fused engine's
#' zero-count decomposition. It runs only for
#' `optop_index_se(reopt = "se", level = "document")`.
#'
#' @keywords internal
.optop_index_se_reopt <- function(theta, phi, dtm, partition, pi_row,
                                  K, macro, ztest, min_null = 0) {
  J <- nrow(theta)
  W <- ncol(phi)
  block_size_doc <- max(100L, min(J, floor(5e8 / (W * 8))))

  D_K <- numeric(J)
  D_null <- numeric(J)
  # deprecated small-corpus path: reconstruct the dense rare rows per block
  # from the sparse partition (the only remaining dense-mask consumer)
  partition <- .optop_partition_upgrade(partition)
  rare_mask <- matrix(FALSE, J, W)
  for (j in seq_len(J)) {
    rare_mask[j, ] <- .optop_partition_rare_row(partition, j, W)
  }

  for (start in seq(1L, J, by = block_size_doc)) {
    end <- min(start + block_size_doc - 1L, J)
    j_idx <- start:end
    block_len <- length(j_idx)

    N_block <- as.matrix(dtm[j_idx, , drop = FALSE])           # block × W
    L_block <- partition$L[j_idx]                              # block vector
    rare_block <- rare_mask[j_idx, , drop = FALSE]             # block × W
    E_block <- (theta[j_idx, , drop = FALSE] %*% phi) * L_block  # block × W
    B_block <- outer(L_block, pi_row)                          # block × W

    # Re-optimization requires per-document lambda computation
    for (local_idx in seq_len(block_len)) {
      j <- j_idx[local_idx]
      rare_j <- rare_block[local_idx, ]

      N_j <- N_block[local_idx, ]
      E_j <- E_block[local_idx, ]
      B_j <- B_block[local_idx, ]

      # Compute lambda for blending E towards B
      Nj_nonrare <- N_j[!rare_j]; Nj_min <- sum(N_j[rare_j])
      Ej_nonrare <- E_j[!rare_j]; Ej_min <- sum(E_j[rare_j])
      Bj_nonrare <- B_j[!rare_j]; Bj_min <- sum(B_j[rare_j])

      num <- sum((Nj_nonrare - Bj_nonrare) * (Ej_nonrare - Bj_nonrare)) +
             (Nj_min - Bj_min) * (Ej_min - Bj_min)
      den <- sum((Ej_nonrare - Bj_nonrare)^2) + (Ej_min - Bj_min)^2
      lambda <- if (den <= 0) 0 else max(0, min(1, num / den))

      # Blend E towards B
      Ej_opt_nonrare <- lambda * Ej_nonrare + (1 - lambda) * Bj_nonrare
      Ej_opt_min <- lambda * Ej_min + (1 - lambda) * Bj_min

      # Compute discrepancies on harmonized support
      sse_k <- sum((Nj_nonrare - Ej_opt_nonrare)^2) + (Nj_min - Ej_opt_min)^2
      sst_se <- sum((Nj_nonrare - Bj_nonrare)^2) + (Nj_min - Bj_min)^2

      D_K[j] <- sse_k
      D_null[j] <- sst_se
    }
  }

  .optop_index_result_doc(D_K, D_null, K, "se", macro, ztest, min_null)
}

#' @describeIn optop_index Pearson chi-square index \eqn{R^2_{chisq}}{R2_chisq},
#'   subject to the min-bin inclusion rule of [`optop_make_partition()`].
#'
#' @export
optop_index_chisq <- function(model, dtm, partition, baseline,
                              macro = FALSE, reopt = c("none", "pearson"),
                              add_baseline_topic = FALSE,
                              level = c("document", "word"),
                              block_size = NULL,
                              ztest = FALSE,
                              n_threads = 1L,
                              min_null = NULL) {

  level <- match.arg(level)
  reopt <- match.arg(reopt)
  min_null <- .optop_resolve_min_null(min_null, partition$c)
  .optop_deprecate_index_args("optop_index_chisq", reopt, add_baseline_topic,
                              ztest)
  if (level == "document") .optop_report_chisq_min(partition)
  res <- .optop_index_chisq_impl(model, dtm, partition, baseline, macro,
                                 reopt, add_baseline_topic, level,
                                 block_size, ztest, n_threads = n_threads,
                                 min_null = min_null)
  if (level == "document") {
    .optop_report_null_floor(res$n_null_excluded, res$null_excluded_share,
                             min_null, "chisq")
  }
  res
}

#' Worker behind optop_index_chisq()
#'
#' See `.optop_index_se_impl()` for the role of these workers and the
#' meaning of `null_disc`.
#'
#' @inheritParams .optop_index_se_impl
#'
#' @return The same list the exported wrapper documents.
#'
#' @keywords internal
.optop_index_chisq_impl <- function(model, dtm, partition, baseline, macro,
                                    reopt, add_baseline_topic, level,
                                    block_size, ztest, null_disc = NULL,
                                    n_threads = 1L, min_null = 0) {

  tp <- optop_as_theta_phi(.optop_materialize_model(model))
  vocab_model <- colnames(tp$phi)
  pi_row <- .optop_validate_alignment(vocab_model, dtm, partition, baseline)

  theta <- tp$theta
  phi <- tp$phi
  J <- nrow(theta)

  if (add_baseline_topic) {
    phi <- rbind(phi, pi_row)
    theta <- cbind(theta, rep(0, J))
  }

  eng <- .optop_index_engine(theta, phi, dtm, partition, pi_row, "chisq",
                             level, block_size, n_threads,
                             do_null = is.null(null_disc))
  if (level == "word") {
    d_null <- if (is.null(null_disc)) eng$chisq$null else null_disc
    return(.optop_index_result_word(eng$chisq$model, d_null, vocab_model,
                                    tp$K, "chisq"))
  }
  D_null <- if (is.null(null_disc)) eng$chisq$null else null_disc
  .optop_index_result_doc(eng$chisq$model, D_null, tp$K, "chisq", macro,
                          ztest, min_null)
}

#' @describeIn optop_index Deviance index \eqn{R^2_{dev}}{R2_dev}, the
#'   default likelihood-based summary: the grouped multinomial deviance at
#'   document level and the Poisson unit deviance at word level. Bounded
#'   above by one; under exact maximum likelihood and nesting the in-sample
#'   document-level index is also nonnegative up to binning effects, a
#'   guarantee that does not extend to approximate inference.
#' @export
optop_index_deviance <- function(model, dtm, partition, baseline,
                                 macro = FALSE, reopt = c("none", "deviance"),
                                 level = c("document", "word"),
                                 block_size = NULL,
                                 ztest = FALSE,
                                 n_threads = 1L,
                                 min_null = NULL) {

  level <- match.arg(level)
  reopt <- match.arg(reopt)
  min_null <- .optop_resolve_min_null(min_null, partition$c)
  .optop_deprecate_index_args("optop_index_deviance", reopt,
                              add_baseline_topic = FALSE, ztest)
  res <- .optop_index_deviance_impl(model, dtm, partition, baseline, macro,
                                    reopt, level, block_size, ztest,
                                    n_threads = n_threads,
                                    min_null = min_null)
  if (level == "document") {
    .optop_report_null_floor(res$n_null_excluded, res$null_excluded_share,
                             min_null, "deviance")
  }
  res
}

#' Worker behind optop_index_deviance()
#'
#' See `.optop_index_se_impl()` for the role of these workers and the
#' meaning of `null_disc`.
#'
#' @inheritParams .optop_index_se_impl
#'
#' @return The same list the exported wrapper documents.
#'
#' @keywords internal
.optop_index_deviance_impl <- function(model, dtm, partition, baseline, macro,
                                       reopt, level, block_size, ztest,
                                       null_disc = NULL, n_threads = 1L,
                                       min_null = 0) {

  tp <- optop_as_theta_phi(.optop_materialize_model(model))
  vocab_model <- colnames(tp$phi)
  pi_row <- .optop_validate_alignment(vocab_model, dtm, partition, baseline)

  theta <- tp$theta
  phi <- tp$phi

  eng <- .optop_index_engine(theta, phi, dtm, partition, pi_row, "deviance",
                             level, block_size, n_threads,
                             do_null = is.null(null_disc))
  if (level == "word") {
    d_null <- if (is.null(null_disc)) eng$deviance$null else null_disc
    return(.optop_index_result_word(eng$deviance$model, d_null, vocab_model,
                                    tp$K, "deviance"))
  }
  D_null <- if (is.null(null_disc)) eng$deviance$null else null_disc
  .optop_index_result_doc(eng$deviance$model, D_null, tp$K, "deviance",
                          macro, ztest, min_null)
}


#' OpTop goodness-of-fit indices table across a grid of LDA models
#'
#' Computes one or more OpTop goodness-of-fit indices as computed by [`optop_index()`],
#' for a list of fitted topic models and returns a tidy table by number of topics \eqn{K}.
#'
#' @param models A non-empty `list` of fitted `topicmodels::LDA` models (VEM or Gibbs),
#'   typically a grid over different \eqn{K}; elements may be loader
#'   functions returning a fit, materialized one at a time.
#' @param dtm A counts document–term matrix (documents × terms). Use
#'   `Matrix::dgCMatrix` or a `quanteda::dfm` that represents raw counts
#'   (not weighted), or an [optop_corpus()] of count shards.
#' @param metrics Character vector selecting which indices to compute. Any subset of
#'   `c("se", "chisq", "deviance")`. Default: `c("se","chisq","deviance")`.
#' @param c Positive scalar used inside [`optop_make_partition()`] to set the per-document
#'   rare threshold \eqn{\tau_j = c / L_j}. Default: `1`, the value the paper
#'   adopts and recommends when the deviance index is primary; use `c = 5`
#'   for Pearson-primary work and report sensitivity to `c`.
#' @param macro Logical; if `TRUE`, also compute macro indices (equal-weight averages of
#'   per-document indices). Default: `FALSE`.
#' @param reopt `r lifecycle::badge("deprecated")` Re-optimization mode
#'   forwarded to the underlying index; scheduled for removal before v1.0.0.
#'   Use the default `"none"`.
#' @param add_baseline_topic `r lifecycle::badge("deprecated")` Logical;
#'   forwarded to the SE and chi-square indices to force non-negativity.
#'   Scheduled for removal before v1.0.0. Default: `FALSE`.
#' @param partition Optional precomputed result from [`optop_make_partition()`].
#'   If `NULL`, it is computed internally from `models` and `dtm`.
#' @param baseline Optional precomputed result from [`optop_make_baseline()`]. If `NULL`,
#'   it is computed internally from `dtm`.
#' @param level Character; aggregation level for the indices. `"document"` (default)
#'   returns document-level micro/macro indices. `"word"` returns word-level
#'   Micro-Word and Macro-Word indices.
#' @param block_size Integer or `NULL`; number of vocabulary terms to process at once
#'   during word-level computation. Smaller values use less memory but may be slower.
#'   If `NULL` (default), automatically chosen based on corpus size to target ~1.5 GB
#'   memory usage. Only used when `level = "word"`.
#' @param ztest `r lifecycle::badge("deprecated")` Logical; if `TRUE`, append
#'   the in-sample Z-test statistics. The methodology reserves inference for
#'   held-out evaluation; scheduled for removal before v1.0.0. Default: `FALSE`.
#' @param n_threads Integer; number of OpenMP threads used by the compiled
#'   index engine and, when `partition` is computed internally, by the
#'   partition kernel (default `1L`). Results are identical for any value.
#' @param min_null Nonnegative scalar; the null-discrepancy floor of the
#'   document-level aggregation, applied uniformly across metrics and models
#'   (the baseline discrepancy does not depend on \eqn{K}, so the same
#'   documents are excluded at every grid point). The default `NULL`
#'   resolves to the partition threshold constant `c`; `min_null = 0`
#'   restores the pre-0.14.1 strict positivity rule. See the
#'   *Null-discrepancy floor* paragraph of [optop_index].
#'
#' @details
#' The function wraps [`optop_index_se()`], [`optop_index_chisq()`], and
#' [`optop_index_deviance()`]. It reuses a harmonized support (via
#' [`optop_make_partition()`]) and a fixed baseline (via
#' [`optop_make_baseline()`]) so results are directly comparable across \eqn{K}.
#'
#' @return When `level = "document"`, a `data.frame` with one row per model and columns:
#' - `K`: number of topics in the model.
#' - `R2_SE`, `R2_chisq`, `R2_dev`: micro indices for the selected metrics.
#' - `R2_SE_macro`, `R2_chisq_macro`, `R2_dev_macro`: macro indices when `macro = TRUE`.
#' - `Z_SE`, `Z_chisq`, `Z_dev`, `pval_SE`, `pval_chisq`, `pval_dev`: Z-test
#'   statistics and p-values when `ztest = TRUE`.
#'
#' When `level = "word"`, a `data.frame` with columns:
#' - `K`: number of topics in the model.
#' - `R2_SE_micro_word`, `R2_chisq_micro_word`, `R2_dev_micro_word`: Micro-Word indices.
#' - `R2_SE_macro_word`, `R2_chisq_macro_word`, `R2_dev_macro_word`: Macro-Word indices.
#'
#' Missing columns are omitted when the corresponding metric is not requested.
#'
#' @section Performance tips:
#' - To avoid recomputation, precompute `partition <- optop_make_partition(...)`
#'   and `baseline <- optop_make_baseline(...)` once and pass them in.
#' - Large vocabularies benefit from the `block` argument inside
#'   [`optop_make_partition()`] (set there, not here).
#'
#' @examples
#' \dontrun{
#' library(quanteda)
#' library(OpTop)
#' 
#' K_grid = 2:10
#' 
#' # Tokenize the corpus
#' toks <- data_corpus_inaugural %>% 
#'   tokens(remove_punct = TRUE, 
#'          remove_symbols = TRUE, 
#'          remove_numbers = TRUE) %>% 
#'   tokens_tolower() %>% 
#'   tokens_remove(stopwords())
#'   
#' # Create the document-feature-matrix
#' mydfm <- dfm(toks)
#' 
#' # Estimate topic models via VEM
#' VEM_models <- lapply(
#'   K_grid, function(k) {
#'     topicmodels::LDA(x = mydfm, k = k)
#'   }
#' )
#'
#' # 1) Quick run: compute partition/baseline internally
#' tab <- optop_index_table(
#'   models   = VEM_models,
#'   dtm      = mydfm,
#'   metrics  = c("se","deviance"),
#'   macro    = TRUE
#' )
#'
#' # 2) Faster repeated runs: precompute and reuse partition/baseline
#' part <- optop_make_partition(models = VEM_models, dtm = mydfm, c = 1)
#' base <- optop_make_baseline(dtm = mydfm)
#' tab2 <- optop_index_table(
#'   models    = VEM_models,
#'   dtm       = mydfm,
#'   metrics   = c("se","chisq","deviance"),
#'   macro     = FALSE,
#'   partition = part,
#'   baseline  = base
#' )
#' }
#'
#' @seealso
#' [`optop_index()`], [`optop_make_partition()`], [`optop_make_baseline()`]
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2026). Goodness-of-fit indices and
#' diagnostics for topic models. Working paper.
#'
#' Lewis, C. M. and Grossetti, F. (2022). A statistical approach for optimal
#' topic model identification. *Journal of Machine Learning Research*,
#' 23(58), 1--20. <https://jmlr.org/papers/v23/19-297.html>
#'
#' @export
optop_index_table <- function(models, dtm, metrics = c("se","chisq","deviance"),
                              c = 1, macro = FALSE, reopt = "none",
                              add_baseline_topic = FALSE,
                              partition = NULL, baseline = NULL,
                              level = c("document", "word"),
                              block_size = NULL,
                              ztest = FALSE,
                              n_threads = 1L,
                              min_null = NULL) {

  level <- match.arg(level)
  stopifnot(length(models) >= 1)
  .optop_deprecate_index_args("optop_index_table", reopt, add_baseline_topic,
                              ztest)
  if (is.null(partition)) {
    partition <- optop_make_partition(models, dtm, c = c,
                                      n_threads = n_threads)
  }
  if (is.null(baseline))  baseline  <- optop_make_baseline(dtm)
  min_null <- .optop_resolve_min_null(min_null, partition$c)

  metrics <- intersect(c("se", "chisq", "deviance"), metrics)
  if ("chisq" %in% metrics && level == "document") {
    .optop_report_chisq_min(partition)
  }
  reopt_se <- if (reopt == "se") "se" else "none"
  # metrics evaluated by the fused engine in one sweep per model; the SE
  # re-optimization path needs the full baseline counts per document and
  # keeps its dedicated implementation (document level only: word-level
  # aggregation never re-optimizes)
  metrics_engine <- if (reopt_se == "se" && level == "document") {
    setdiff(metrics, "se")
  } else {
    metrics
  }

  # The baseline (no-topics) discrepancy does not depend on the fitted
  # model: one engine sweep computes it for every metric, shared across the
  # whole grid of K.
  dtm_vocab <- if (.optop_is_corpus(dtm)) {
    .optop_corpus_vocab(dtm)
  } else {
    colnames(dtm)
  }
  pi_row_tab <- .optop_validate_alignment(dtm_vocab, dtm, partition,
                                          baseline)
  nulls <- if (length(metrics_engine)) {
    .optop_index_engine(theta = NULL, phi = NULL, dtm, partition, pi_row_tab,
                        metrics_engine, level, block_size, n_threads,
                        do_model = FALSE, do_null = TRUE)
  } else {
    NULL
  }

  rows <- vector("list", length(models))
  # D_null is model-independent, so the floor excludes the same documents at
  # every K: one report per metric, taken from the first model's result
  floor_report <- list()
  for (i_mod in seq_along(models)) {
    m <- .optop_materialize_model(models[[i_mod]])
    tp <- optop_as_theta_phi(m)
    vocab_model <- colnames(tp$phi)
    pi_row <- .optop_validate_alignment(vocab_model, dtm, partition, baseline)
    res <- list(K = tp$K)

    # one fused sweep per model for every engine-evaluated metric; the
    # baseline augmentation is inert under reopt = "none" (a zero-weight
    # topic leaves the fitted counts unchanged), so one set of matrices
    # serves all metrics
    theta <- tp$theta
    phi <- tp$phi
    if (add_baseline_topic) {
      phi <- rbind(phi, pi_row)
      theta <- cbind(theta, rep(0, nrow(theta)))
    }
    eng <- if (length(metrics_engine)) {
      .optop_index_engine(theta, phi, dtm, partition, pi_row,
                          metrics_engine, level, block_size, n_threads,
                          do_model = TRUE, do_null = FALSE)
    } else {
      NULL
    }

    if (level == "word") {
      if ("se" %in% metrics) {
        x <- .optop_index_result_word(eng$se$model, nulls$se$null,
                                      vocab_model, tp$K, "se")
        res$R2_SE_micro_word <- x$r2_micro_word
        res$R2_SE_macro_word <- x$r2_macro_word
      }
      if ("chisq" %in% metrics) {
        x <- .optop_index_result_word(eng$chisq$model, nulls$chisq$null,
                                      vocab_model, tp$K, "chisq")
        res$R2_chisq_micro_word <- x$r2_micro_word
        res$R2_chisq_macro_word <- x$r2_macro_word
      }
      if ("deviance" %in% metrics) {
        x <- .optop_index_result_word(eng$deviance$model,
                                      nulls$deviance$null,
                                      vocab_model, tp$K, "deviance")
        res$R2_dev_micro_word <- x$r2_micro_word
        res$R2_dev_macro_word <- x$r2_macro_word
      }
    } else {
      if ("se" %in% metrics) {
        x <- if (reopt_se == "se") {
          .optop_index_se_impl(m, dtm, partition, baseline,
                               macro = macro, reopt = "se",
                               add_baseline_topic = add_baseline_topic,
                               level = "document", block_size = block_size,
                               ztest = ztest, n_threads = n_threads,
                               min_null = min_null)
        } else {
          .optop_index_result_doc(eng$se$model, nulls$se$null, tp$K, "se",
                                  macro, ztest, min_null)
        }
        if (i_mod == 1L) floor_report$se <- x
        res$R2_SE <- x$r2
        if (macro) res$R2_SE_macro <- x$r2_macro
        if (ztest && !is.null(x$ztest)) {
          res$Z_SE <- x$ztest$z
          res$pval_SE <- x$ztest$pval
        }
      }
      if ("chisq" %in% metrics) {
        x <- .optop_index_result_doc(eng$chisq$model, nulls$chisq$null,
                                     tp$K, "chisq", macro, ztest, min_null)
        if (i_mod == 1L) floor_report$chisq <- x
        res$R2_chisq <- x$r2
        if (macro) res$R2_chisq_macro <- x$r2_macro
        if (ztest && !is.null(x$ztest)) {
          res$Z_chisq <- x$ztest$z
          res$pval_chisq <- x$ztest$pval
        }
      }
      if ("deviance" %in% metrics) {
        x <- .optop_index_result_doc(eng$deviance$model, nulls$deviance$null,
                                     tp$K, "deviance", macro, ztest,
                                     min_null)
        if (i_mod == 1L) floor_report$deviance <- x
        res$R2_dev <- x$r2
        if (macro) res$R2_dev_macro <- x$r2_macro
        if (ztest && !is.null(x$ztest)) {
          res$Z_dev <- x$ztest$z
          res$pval_dev <- x$ztest$pval
        }
      }
    }
    rows[[i_mod]] <- as.data.frame(res, check.names = FALSE)
  }
  for (metric in names(floor_report)) {
    x <- floor_report[[metric]]
    .optop_report_null_floor(x$n_null_excluded, x$null_excluded_share,
                             min_null, metric)
  }
  out <- do.call(rbind, rows)
  data.table::setDT(out)[]
}


#' Fused compiled evaluation of the discrepancy indices
#'
#' Routes both aggregation levels of every index through the compiled
#' kernels of `src/index_core.cpp`, which reproduce the conventions of the
#' original block-vectorized R implementation exactly (per-metric eps
#' floors, collapsed-bin sums, decomposition into full, rare and min terms)
#' while computing all requested metrics, and optionally the
#' model-independent baseline side, in a single traversal with no
#' `J x W` temporaries. Setting `do_model = FALSE` evaluates the baseline
#' side alone, which is how the grid-level `D_null` hoist is computed.
#'
#' @param theta,phi The (possibly baseline-augmented) model matrices;
#'   ignored when `do_model = FALSE`.
#' @param dtm The counts document-term matrix.
#' @param partition A partition from [optop_make_partition()].
#' @param pi_row Baseline probabilities aligned to the dtm vocabulary.
#' @param metrics Character subset of `c("se", "chisq", "deviance")`.
#' @param level Aggregation level, `"document"` or `"word"`.
#' @param block_size Words per block at the word level (`NULL` for the
#'   adaptive default); the document level derives its own block size.
#' @param n_threads OpenMP threads; never changes the result.
#' @param do_model,do_null Which sides to evaluate.
#'
#' @return A list with one entry per metric (`se`, `chisq`, `deviance`),
#'   each a list of numeric vectors `model` and `null` (per document at the
#'   document level, per word at the word level); metrics or sides that were
#'   not requested are zero vectors.
#'
#' @keywords internal
.optop_index_engine <- function(theta, phi, dtm, partition, pi_row,
                                metrics, level, block_size = NULL,
                                n_threads = 1L, do_model = TRUE,
                                do_null = TRUE) {
  # sharded corpora stream one shard at a time: per-document rows (document
  # level) and per-word accumulators (word level) combine exactly across
  # shards because every document is independent
  if (.optop_is_corpus(dtm)) {
    # word level needs only the document lengths from the partition object,
    # so length-only stubs (the fit-strata instrument builder) stay valid
    if (level == "document") {
      partition <- .optop_partition_upgrade(partition)
    }
    acc <- NULL
    offset <- 0
    for (s in seq_len(dtm$n_shards)) {
      m <- .optop_corpus_shard(dtm, s)
      n_s <- nrow(m)
      idx <- offset + seq_len(n_s)
      part_s <- if (level == "document") {
        .optop_partition_slice(partition, offset, n_s)
      } else {
        list(L = partition$L[idx])
      }
      theta_s <- if (do_model) theta[idx, , drop = FALSE] else theta
      eng_s <- .optop_index_engine(theta_s, phi, m, part_s, pi_row, metrics,
                                   level, block_size, n_threads,
                                   do_model, do_null)
      if (level == "word") {
        if (is.null(acc)) {
          acc <- eng_s
        } else {
          for (metric in names(acc)) {
            acc[[metric]]$model <- acc[[metric]]$model + eng_s[[metric]]$model
            acc[[metric]]$null <- acc[[metric]]$null + eng_s[[metric]]$null
          }
        }
      } else {
        if (is.null(acc)) {
          J_tot <- length(partition$L)
          acc <- lapply(eng_s, function(e) {
            list(model = numeric(J_tot), null = numeric(J_tot))
          })
        }
        for (metric in names(acc)) {
          acc[[metric]]$model[idx] <- eng_s[[metric]]$model
          acc[[metric]]$null[idx] <- eng_s[[metric]]$null
        }
      }
      offset <- offset + n_s
    }
    return(acc)
  }

  do_se <- "se" %in% metrics
  do_chisq <- "chisq" %in% metrics
  do_dev <- "deviance" %in% metrics
  eps <- 1e-12
  N <- methods::as(dtm, "CsparseMatrix")
  J <- nrow(N)
  W <- ncol(N)
  L <- as.numeric(partition$L)
  pi_num <- as.numeric(pi_row)
  n_threads <- as.integer(n_threads)

  if (level == "word") {
    # the document dimension is blocked INSIDE the kernel (a block x w_len
    # gemm buffer), so the word block size only bounds the phi slice and the
    # accumulator width; counts cross the boundary as zero-copy CSC slots
    if (is.null(block_size)) {
      block_size <- 4096L
    }
    theta_in <- if (do_model) theta else matrix(0, 0, 0)
    acc <- matrix(0, W, 6)
    for (start in seq(1L, W, by = block_size)) {
      end <- min(start + block_size - 1L, W)
      w_idx <- start:end
      phi_cols <- if (do_model) {
        phi[, w_idx, drop = FALSE]
      } else {
        matrix(0, 0, 0)
      }
      acc[w_idx, ] <- optop_index_word_core(theta_in, phi_cols,
                                            N@p, N@i, N@x,
                                            start - 1L, length(w_idx), L,
                                            pi_num[w_idx], eps,
                                            do_model, do_null,
                                            do_se, do_chisq, do_dev,
                                            n_threads)
    }
  } else {
    # sparse-partition merge-join: one call over all documents, each at
    # O(|NR_j| * K + nnz_j); dense-mask partitions are upgraded on entry
    partition <- .optop_partition_upgrade(partition)
    N_t <- Matrix::t(N)
    # Pearson min-bin inclusion flags: partitions from OpTop < 0.13.0 lack
    # the field and keep the pre-rule behavior (every min bin included)
    min_ok <- partition$chisq_min_ok
    if (is.null(min_ok)) {
      if (do_chisq) {
        cli::cli_warn(paste(
          "the partition has no {.field chisq_min_ok} field (created by",
          "OpTop < 0.13.0); the Pearson min bin is kept for every document.",
          "Recompute {.fun optop_make_partition} to apply the inclusion rule."
        ))
      }
      min_ok <- rep(TRUE, J)
    }
    theta_in <- if (do_model) theta else matrix(0, 0, 0)
    phi_in <- if (do_model) phi else matrix(0, 0, 0)
    acc <- optop_index_doc_core(theta_in, phi_in, N_t@p, N_t@i, N_t@x,
                                partition$nonrare_offsets,
                                partition$nonrare_words,
                                L, pi_num, min_ok, eps,
                                do_model, do_null,
                                do_se, do_chisq, do_dev,
                                n_threads)
  }
  list(se = list(model = acc[, 1], null = acc[, 2]),
       chisq = list(model = acc[, 3], null = acc[, 4]),
       deviance = list(model = acc[, 5], null = acc[, 6]))
}

# Assemble the document-level result list from the per-document
# discrepancies. Aggregation is restricted to the non-degenerate documents
# J+ = {j : D_j(null) >= min_null}: the Micro sums run over J+ only and
# excluded documents carry an undefined (NA) document-level index. With
# min_null = 0 the set falls back to the legacy strict positivity rule
# J+ = {j : D_j(null) > 0}, which reproduces pre-0.14.1 output exactly.
# The floor is the family-generic regularity condition on the Macro
# estimand: at partition resolution tau_j = c/L_j a document whose support
# is swept into the min bin has observed == baseline by construction, so
# D_j(null) ~ 0 in every family at once and r_j = 1 - D_j(K)/D_j(null) is
# unbounded below. Micro is nearly immune (weights prop. to D_j(null)).
.optop_index_result_doc <- function(D_K, D_null, K, metric, macro, ztest,
                                    min_null = 0) {
  valid <- if (min_null > 0) D_null >= min_null else D_null > 0
  r2_doc <- rep(NA_real_, length(D_K))
  r2_doc[valid] <- 1 - D_K[valid] / D_null[valid]
  r2_micro <- 1 - sum(D_K[valid]) / sum(D_null[valid])
  r2_macro <- mean(r2_doc[valid])
  result <- list(r2 = r2_micro, r2_macro = if (macro) r2_macro else NULL,
                 r2_doc = r2_doc, d_model = D_K, d_null = D_null,
                 n_null_excluded = sum(!valid & is.finite(D_null)),
                 null_excluded_share =
                   sum(!valid & is.finite(D_null)) / length(D_null),
                 K = K, metric = metric)
  if (ztest) {
    result$ztest <- .optop_ztest(r2_doc[valid], r2_macro)
  }
  result
}

# Assemble the word-level result list from the per-word discrepancies. Words
# with zero baseline discrepancy have no defined index (NA); both summaries
# aggregate over V+ = {w : d_null(w) > 0} only.
.optop_index_result_word <- function(d_w, d_null, vocab, K, metric) {
  valid <- d_null > 0
  r2_word <- rep(NA_real_, length(d_w))
  r2_word[valid] <- 1 - d_w[valid] / d_null[valid]
  names(r2_word) <- vocab
  omega <- d_null[valid] / sum(d_null[valid])
  list(r2_word = r2_word,
       r2_micro_word = sum(omega * r2_word[valid]),
       r2_macro_word = mean(r2_word[valid]),
       d_model = d_w, d_null = d_null,
       K = K, metric = metric)
}

#' Validate the alignment of dtm, partition and baseline with a model
#'
#' Shared gate of the index functions: the dtm, the baseline and the
#' harmonized partition must all match the model vocabulary (same features,
#' same order); mismatches stop with a pointer to
#' `optop_align_dtm_to_models()`.
#'
#' @param vocab_model Character vector; the model's vocabulary in model
#'   order.
#' @param dtm The counts document-term matrix under validation.
#' @param partition A partition from [optop_make_partition()].
#' @param baseline A baseline from [optop_make_baseline()].
#'
#' @return The baseline probabilities ordered as `vocab_model`.
#'
#' @keywords internal
.optop_validate_alignment <- function(vocab_model, dtm, partition, baseline) {
  dtm_vocab <- if (.optop_is_corpus(dtm)) {
    .optop_corpus_vocab(dtm)
  } else {
    colnames(dtm)
  }
  if (!identical(dtm_vocab, vocab_model))
    stop("DTM vocabulary/order differs from model. Use optop_align_dtm_to_models() and recompute partition/baseline.")
  # baseline alignment
  if (!is.null(names(baseline$pi_glob))) {
    pi_row <- baseline$pi_glob[vocab_model]
    if (any(is.na(pi_row)))
      stop("Baseline vocabulary does not match model. Recompute optop_make_baseline() on an aligned DTM.")
  } else {
    if (length(baseline$pi_glob) != length(vocab_model))
      stop("Baseline length != model vocab. Recompute baseline on aligned DTM.")
    pi_row <- baseline$pi_glob
  }

  # partition alignment: sparse partitions carry the vocabulary explicitly,
  # dense-mask partitions (OpTop < 0.15.0) carry it as mask column names
  part_vocab <- if (identical(partition$format, 2L)) {
    partition$vocab
  } else {
    colnames(partition$rare_mask)
  }
  if (is.null(part_vocab) || !identical(part_vocab, vocab_model)) {
    stop("Partition vocabulary != model vocabulary. Recompute optop_make_partition() on the aligned DTM.")
  }
  pi_row
}

#' Baseline (no-topics) discrepancy for one metric
#'
#' This quantity does not depend on the fitted model, only on the observed
#' counts, the harmonized partition and the global baseline, so grid
#' evaluations across \eqn{K} can compute it once and pass it to the index
#' workers via their `null_disc` argument. Formulas, eps floors and blocking
#' mirror the corresponding index workers exactly.
#'
#' @param dtm The counts document-term matrix.
#' @param partition A partition from [optop_make_partition()].
#' @param pi_row Baseline probabilities aligned to the dtm vocabulary, as
#'   returned by `.optop_validate_alignment()`.
#' @param metric The discrepancy metric.
#' @param level Aggregation level, `"document"` or `"word"`.
#' @param block_size Number of vocabulary columns per processing block.
#' @param n_threads OpenMP threads used by the fused engine.
#'
#' @return The per-document \eqn{D_{null}}{D_null} vector when
#'   `level = "document"` and the per-word null discrepancy vector when
#'   `level = "word"`.
#'
#' @keywords internal
.optop_index_null <- function(dtm, partition, pi_row,
                              metric = c("se", "chisq", "deviance"),
                              level = c("document", "word"),
                              block_size = NULL, n_threads = 1L) {
  metric <- match.arg(metric)
  level <- match.arg(level)
  eng <- .optop_index_engine(theta = NULL, phi = NULL, dtm, partition,
                             pi_row, metric, level, block_size, n_threads,
                             do_model = FALSE, do_null = TRUE)
  eng[[metric]]$null
}

#' Deprecation gate for the non-paper index arguments
#'
#' Emits the lifecycle warnings for the three features that Lewis and
#' Grossetti (2026) does not support: the re-optimization blending, the
#' baseline-topic augmentation (both enforce non-negativity, whereas the
#' methodology treats negative indices as informative), and the in-sample
#' Z-test (the methodology reserves inference for held-out evaluation).
#' Behavior is unchanged until removal before v1.0.0.
#'
#' @param fn Name of the calling exported function, for the warning text.
#' @param reopt,add_baseline_topic,ztest The argument values to inspect.
#'
#' @return Invisibly `NULL`; called for its side effects.
#'
#' @keywords internal
.optop_deprecate_index_args <- function(fn, reopt, add_baseline_topic,
                                        ztest) {
  if (!identical(reopt, "none")) {
    lifecycle::deprecate_warn(
      when = "0.13.0",
      what = paste0(fn, "(reopt)"),
      details = paste("The re-optimization blending is not part of the",
                      "methodology of Lewis and Grossetti (2026) and will be",
                      "removed before v1.0.0.")
    )
  }
  if (isTRUE(add_baseline_topic)) {
    lifecycle::deprecate_warn(
      when = "0.13.0",
      what = paste0(fn, "(add_baseline_topic)"),
      details = paste("The augmentation enforces non-negativity, while the",
                      "methodology treats negative indices as informative;",
                      "it will be removed before v1.0.0.")
    )
  }
  if (isTRUE(ztest)) {
    lifecycle::deprecate_warn(
      when = "0.13.0",
      what = paste0(fn, "(ztest)"),
      details = paste("The methodology reserves inference for held-out",
                      "evaluation; the in-sample Z-test will be removed",
                      "before v1.0.0.")
    )
  }
  invisible(NULL)
}

#' Resolve the null-discrepancy floor of the document-level indices
#'
#' The default `min_null = NULL` resolves to the partition threshold
#' constant `c`, so the conditioning event of the document-level
#' aggregation sharpens from \eqn{\{D_j(\mathrm{null}) > 0\}}{{D_j(null) > 0}}
#' to \eqn{\{D_j(\mathrm{null}) \ge c\}}{{D_j(null) >= c}}. Partitions
#' created before the field existed resolve to `0`, the legacy rule.
#'
#' @param min_null The user-supplied value (`NULL` for the default).
#' @param c_part The partition's threshold constant (`partition$c`).
#'
#' @return A single nonnegative number.
#'
#' @keywords internal
.optop_resolve_min_null <- function(min_null, c_part) {
  if (is.null(min_null)) {
    return(if (is.null(c_part)) 0 else c_part)
  }
  if (!is.numeric(min_null) || length(min_null) != 1L ||
      !is.finite(min_null) || min_null < 0) {
    stop("min_null must be a single nonnegative number")
  }
  min_null
}

#' Report the null-discrepancy floor exclusions
#'
#' Signals once per call how many documents were excluded from the
#' document-level aggregation because \eqn{D_j(\mathrm{null}) <
#' \mathrm{min\_null}}{D_j(null) < min_null}. Tiny baseline discrepancies
#' signal structurally collapsed supports rather than sampling noise: under
#' the null a document's binned baseline discrepancy is approximately
#' chi-square with (number of effective bins - 1) degrees of freedom, so
#' values orders of magnitude below 1 can only arise when the partition has
#' degenerated to essentially one bin.
#'
#' @param n_excluded Number of excluded documents.
#' @param share Excluded share of the corpus.
#' @param min_null The resolved floor.
#' @param metric The discrepancy family, for the message.
#'
#' @return Invisibly `NULL`; called for its side effect.
#'
#' @keywords internal
.optop_report_null_floor <- function(n_excluded, share, min_null, metric) {
  # min_null = 0 is the legacy escape hatch: pre-0.14.1 behavior, silence
  # included (documents with D_null exactly zero were always excluded)
  if (min_null <= 0 || !isTRUE(n_excluded > 0)) return(invisible(NULL))
  cli::cli_alert_info(paste(
    "{n_excluded} document{?s} excluded from the {metric} aggregation by",
    "the null-discrepancy floor (D_null < {min_null};",
    "{sprintf('%.1f%%', 100 * share)} of the corpus)"
  ))
  invisible(NULL)
}

#' Report the Pearson min-bin exclusions of a partition
#'
#' Implements the reporting side of the inclusion rule: when any document's
#' collapsed min-bin fails
#' \eqn{\min(\min_K E^K_{j,\min}, B_{j,\min}) \ge c}{min(min_K E_min, B_min) >= c},
#' the share of documents affected and the mean excluded observed mass are
#' signalled once per call.
#'
#' @param partition A partition from [optop_make_partition()].
#'
#' @return Invisibly `NULL`; called for its side effect.
#'
#' @keywords internal
.optop_report_chisq_min <- function(partition) {
  rep <- partition$chisq_min_report
  if (is.null(rep) || !isTRUE(rep$n_excluded > 0)) return(invisible(NULL))
  cli::cli_alert_info(paste(
    "Pearson min bin excluded for {rep$n_excluded} document{?s}",
    "({sprintf('%.1f%%', 100 * rep$share)} of the corpus; mean excluded",
    "observed mass {sprintf('%.4f', rep$excluded_mass)})"
  ))
  invisible(NULL)
}

#' Z-test for cross-document inference on the Macro R-squared index
#'
#' Tests whether the \eqn{K}-topic model provides statistically significant
#' improvement over the no-topics baseline.
#'
#' @param r2_doc Numeric vector of document-level \eqn{R^2}{R^2} values
#'   (length \eqn{J}).
#' @param r2_macro Scalar macro index (mean of `r2_doc` over valid
#'   documents).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{z}: Z-statistic.
#'   \item \code{pval}: One-sided p-value for
#'     \eqn{H_1: \mu_{R^2} > 0}{H1: mu_R2 > 0}.
#'   \item \code{se}: Standard error \eqn{\hat\sigma_R}{sigma_R hat}.
#'   \item \code{ci}: 95\% confidence interval for the true mean
#'     \eqn{R^2}{R^2}.
#'   \item \code{J}: Number of valid documents used.
#' }
#'
#' @details
#' The null hypothesis is \eqn{H_0: \mu_{R^2} \le 0}{H0: mu_R2 <= 0} (the
#' topic model is no better than the global distribution on average). Under
#' regularity conditions, the test statistic
#' \eqn{Z = \sqrt{J}\,\bar{R}^2_{Macro} / \hat\sigma_R}{Z = sqrt(J) * R2_Macro / sigma_R hat}
#' is asymptotically \eqn{N(0, 1)}.
#'
#' @keywords internal
.optop_ztest <- function(r2_doc, r2_macro) {
  # Exclude degenerate documents (where baseline discrepancy = 0)
  valid <- !is.na(r2_doc) & is.finite(r2_doc)
  r2_valid <- r2_doc[valid]
  J <- length(r2_valid)

  if (J < 2) {
    return(list(z = NA_real_, pval = NA_real_, se = NA_real_,
                ci = c(NA_real_, NA_real_), J = J))
  }

  # Variance estimator: σ̂²_R = (1/(J-1)) Σ_j (R²_j - R̄²_Macro)²
  sigma2_hat <- sum((r2_valid - r2_macro)^2) / (J - 1)
  sigma_hat <- sqrt(sigma2_hat)

  # Z-statistic: Z = √J · R̄²_Macro / σ̂_R
  se_mean <- sigma_hat / sqrt(J)
  z_stat <- r2_macro / se_mean

  # One-sided p-value for H1: μ_R² > 0
  pval <- stats::pnorm(z_stat, lower.tail = FALSE)

  # 95% CI for the true mean
  ci <- r2_macro + c(-1, 1) * stats::qnorm(0.975) * se_mean

  list(z = z_stat, pval = pval, se = se_mean, ci = ci, J = J)
}
