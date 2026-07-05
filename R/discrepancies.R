#' OpTop goodness-of-fit indices for a single LDA model
#'
#' Compute regression-style goodness-of-fit indices based on sum of squared errors (SE),
#' chi-square, and deviance for one fitted topic model. Each index compares
#' model-expected counts to a fixed baseline (the corpus word distribution) on a
#' harmonized support per document, where rare words are collapsed into one bin as
#' defined by [`optop_make_partition()`].
#'
#' @param model A single fitted `topicmodels::LDA` model (VEM or Gibbs).
#' @param dtm A counts document–term matrix (documents × terms). Use `Matrix::dgCMatrix`
#'   or a `quanteda::dfm` that represents raw counts (not weighted).
#' @param partition Result from [`optop_make_partition()`], computed on a DTM that is
#'   *aligned* to the model vocabulary (same features and order).
#' @param baseline Result from [`optop_make_baseline()`], computed on the same aligned
#'   counts DTM used for `partition`.
#' @param macro Logical; if `TRUE`, also compute the macro index (equal-weight average of
#'   the per-document indices). Default: `FALSE`.
#' @param reopt Character scalar indicating optional re-optimization; routed per metric.
#'   Allowed values include `"none"` (default), `"se"`, `"pearson"`, and `"deviance"`.
#'   Unsupported values for a given metric are ignored.
#' @param level Character; aggregation level for the index. `"document"` (default) computes
#'   document-level indices aggregated across words. `"word"` computes word-level indices
#'   aggregated across documents.
#' @param block_size Integer or `NULL`; number of vocabulary terms to process at once
#'   during word-level computation. Smaller values use less memory but may be slower.
#'   If `NULL` (default), automatically chosen based on corpus size to target ~1.5 GB
#'   memory usage. Only used when `level = "word"`.
#' @param ztest Logical; if `TRUE`, append a Z-test for cross-document inference.
#'   Tests whether the topic model provides statistically significant improvement
#'   over the no-topics baseline. Default: `FALSE`. Only applicable when `level = "document"`.
#' @param n_threads Integer; number of OpenMP threads used by the compiled
#'   index engine (default `1L`). Results are identical for any value; only
#'   wall time changes.
#'
#' @details
#' **Harmonized support.** For each document, rare words (as determined by
#' [`optop_make_partition()`]) are collapsed into a single "min" bin. Indices are then
#' evaluated on `{non-rare terms} ∪ {min}` to ensure comparability across \eqn{K}.
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
#' **Document-level aggregation** (`level = "document"`). Returns the micro index
#' (e.g., \eqn{R^2_{SE}}{R2_SE}) and, if requested, the macro index
#' (e.g., \eqn{\bar R^2_{SE}}{R2_SE_bar}), plus per-document components
#' (e.g., \eqn{R^2_{SE,j}}{R2_SE,j}).
#'
#' **Word-level aggregation** (`level = "word"`). Returns per-word indices
#' \eqn{R^2_{D,w}(K)}, plus Micro-Word (frequency-weighted) and Macro-Word (unweighted)
#' corpus-level summaries. This perspective reveals which words are well-captured vs.
#' poorly modeled by the topic structure.
#'
#' **Z-test** (`ztest = TRUE`). Implements a hypothesis test for cross-document inference.
#' Under H0: \eqn{\mu_{R^2} \leq 0} (no improvement), the statistic \eqn{Z = \sqrt{J}\cdot \bar{R}^2_{Macro} / \hat{\sigma}_R}R
#' is asymptotically \eqn{N(0,1)}. Requires `macro = TRUE` implicitly.
#'
#' @return When `level = "document"`, a list with:
#' - `r2`: scalar micro index (e.g., \eqn{R^2_{SE}}{R2_SE}).
#' - `r2_macro`: scalar macro index if `macro = TRUE`, otherwise `NULL`.
#' - `r2_doc`: numeric vector (length J) with per-document contributions.
#' - `K`: number of topics in `model`.
#' - `metric`: one of `"se"`, `"chisq"`, `"deviance"`.
#' - `ztest`: (if `ztest = TRUE`) list with `z`, `pval`, `se`, `ci`, `J`.
#'
#' When `level = "word"`, a list with:
#' - `r2_word`: named numeric vector (length W) of per-word \eqn{R^2_{D,w}(K)}.
#' - `r2_micro_word`: scalar Micro-Word index (frequency-weighted average).
#' - `r2_macro_word`: scalar Macro-Word index (unweighted average).
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
#' partition <- optop_make_partition(models = VEM_models, dtm = mydfm, c = 5)
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
#' @param add_baseline_topic Logical; if `TRUE`, augment topics with a baseline row.
#'   Combined with `reopt = "se"`, this guarantees non-negativity of \eqn{R^2}.
#'   With `reopt = "none"` (default), the augmentation has no effect on results.
#'   (Only for `optop_index_se()` and `optop_index_chisq()`).
#' @export
optop_index_se <- function(model, dtm, partition, baseline,
                           macro = FALSE, reopt = c("none", "se"),
                           add_baseline_topic = TRUE,
                           level = c("document", "word"),
                           block_size = NULL,
                           ztest = FALSE,
                           n_threads = 1L) {

  level <- match.arg(level)
  reopt <- match.arg(reopt)
  .optop_index_se_impl(model, dtm, partition, baseline, macro, reopt,
                       add_baseline_topic, level, block_size, ztest,
                       n_threads = n_threads)
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
                                 null_disc = NULL, n_threads = 1L) {

  tp <- optop_as_theta_phi(model)
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
                                 tp$K, macro, ztest))
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
  .optop_index_result_doc(eng$se$model, D_null, tp$K, "se", macro, ztest)
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
                                  K, macro, ztest) {
  J <- nrow(theta)
  W <- ncol(phi)
  block_size_doc <- max(100L, min(J, floor(5e8 / (W * 8))))

  D_K <- numeric(J)
  D_null <- numeric(J)
  r2_doc <- numeric(J)
  rare_mask <- partition$rare_mask

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
      r2_doc[j] <- if (sst_se > 0) 1 - sse_k / sst_se else 0
    }
  }

  r2_micro <- 1 - sum(D_K) / sum(D_null)
  r2_macro <- mean(r2_doc[D_null > 0])
  result <- list(r2 = r2_micro, r2_macro = if (macro) r2_macro else NULL,
                 r2_doc = r2_doc, K = K, metric = "se")
  if (ztest) {
    result$ztest <- .optop_ztest(r2_doc[D_null > 0], r2_macro)
  }
  result
}

#' @describeIn optop_index Pearson chi-square index \eqn{R^2_{chisq}}{R2_chisq}.
#' @param add_baseline_topic Logical; if `TRUE`, augment topics with a baseline row.
#'   Combined with `reopt = "se"`, this guarantees non-negativity of \eqn{R^2}.
#'   With `reopt = "none"` (default), the augmentation has no effect on results.
#' (Only for `optop_index_se()` and `optop_index_chisq()`).
#'
#' @export
optop_index_chisq <- function(model, dtm, partition, baseline,
                              macro = FALSE, reopt = c("none", "pearson"),
                              add_baseline_topic = TRUE,
                              level = c("document", "word"),
                              block_size = NULL,
                              ztest = FALSE,
                              n_threads = 1L) {

  level <- match.arg(level)
  reopt <- match.arg(reopt)
  .optop_index_chisq_impl(model, dtm, partition, baseline, macro, reopt,
                          add_baseline_topic, level, block_size, ztest,
                          n_threads = n_threads)
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
                                    n_threads = 1L) {

  tp <- optop_as_theta_phi(model)
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
  .optop_index_result_doc(eng$chisq$model, D_null, tp$K, "chisq", macro, ztest)
}

#' @describeIn optop_index Deviance index \eqn{R^2_{dev}}{R2_dev}.
#'
#' Based on the multinomial deviance; monotone in \eqn{K} under (approximate) ML fits.
#' @export
optop_index_deviance <- function(model, dtm, partition, baseline,
                                 macro = FALSE, reopt = c("none", "deviance"),
                                 level = c("document", "word"),
                                 block_size = NULL,
                                 ztest = FALSE,
                                 n_threads = 1L) {

  level <- match.arg(level)
  reopt <- match.arg(reopt)
  .optop_index_deviance_impl(model, dtm, partition, baseline, macro, reopt,
                             level, block_size, ztest,
                             n_threads = n_threads)
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
                                       null_disc = NULL, n_threads = 1L) {

  tp <- optop_as_theta_phi(model)
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
                          macro, ztest)
}


#' OpTop goodness-of-fit indices table across a grid of LDA models
#'
#' Computes one or more OpTop goodness-of-fit indices as computed by [`optop_index()`],
#' for a list of fitted topic models and returns a tidy table by number of topics \eqn{K}.
#'
#' @param models A non-empty `list` of fitted `topicmodels::LDA` models (VEM or Gibbs),
#'   typically a grid over different \eqn{K}.
#' @param dtm A counts document–term matrix (documents × terms). Use
#'   `Matrix::dgCMatrix` or a `quanteda::dfm` that represents raw counts
#'   (not weighted).
#' @param metrics Character vector selecting which indices to compute. Any subset of
#'   `c("se", "chisq", "deviance")`. Default: `c("se","chisq","deviance")`.
#' @param c Positive scalar used inside [`optop_make_partition()`] to set the per-document
#'   rare threshold \eqn{\tau_j = c / L_j}. Default: `5`.
#' @param macro Logical; if `TRUE`, also compute macro indices (equal-weight averages of
#'   per-document indices). Default: `FALSE`.
#' @param reopt Re-optimization mode forwarded to the underlying index:
#'   - `"none"` (default): no re-optimization;
#'   - `"se"`: enable SE blending in [`optop_index_se()`];
#'   - `"deviance"`: reserved for a deviance-specific re-optimizer (if implemented).
#'
#'   Values are routed per metric (e.g., `"se"` only affects the SE column).
#' @param add_baseline_topic Logical; forwarded to SE and chi-square indices to augment
#'   the topics with a baseline row and guarantee non-negativity of their \eqn{R^2}.
#'   Ignored by the deviance index. Default: `TRUE`.
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
#' @param ztest Logical; if `TRUE`, append Z-test statistics for cross-document inference.
#'   Only applicable when `level = "document"`. Default: `FALSE`.
#' @param n_threads Integer; number of OpenMP threads used by the compiled
#'   index engine and, when `partition` is computed internally, by the
#'   partition kernel (default `1L`). Results are identical for any value.
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
#'   macro    = TRUE,
#'   reopt    = "se"          # enable SE blending only
#' )
#'
#' # 2) Faster repeated runs: precompute and reuse partition/baseline
#' part <- optop_make_partition(models = VEM_models, dtm = mydfm, c = 5)
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
#' Lewis, C. M. and Grossetti, F. (2022). A statistical approach for optimal
#' topic model identification. *Journal of Machine Learning Research*,
#' 23(58), 1--20. <https://jmlr.org/papers/v23/19-297.html>
#'
#' @export
optop_index_table <- function(models, dtm, metrics = c("se","chisq","deviance"),
                              c = 5, macro = FALSE, reopt = "none",
                              add_baseline_topic = TRUE,
                              partition = NULL, baseline = NULL,
                              level = c("document", "word"),
                              block_size = NULL,
                              ztest = FALSE,
                              n_threads = 1L) {

  level <- match.arg(level)
  stopifnot(length(models) >= 1)
  if (is.null(partition)) {
    partition <- optop_make_partition(models, dtm, c = c,
                                      n_threads = n_threads)
  }
  if (is.null(baseline))  baseline  <- optop_make_baseline(dtm)

  metrics <- intersect(c("se", "chisq", "deviance"), metrics)
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
  pi_row_tab <- .optop_validate_alignment(colnames(dtm), dtm, partition,
                                          baseline)
  nulls <- if (length(metrics_engine)) {
    .optop_index_engine(theta = NULL, phi = NULL, dtm, partition, pi_row_tab,
                        metrics_engine, level, block_size, n_threads,
                        do_model = FALSE, do_null = TRUE)
  } else {
    NULL
  }

  rows <- vector("list", length(models))
  for (i_mod in seq_along(models)) {
    m <- models[[i_mod]]
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
                               ztest = ztest, n_threads = n_threads)
        } else {
          .optop_index_result_doc(eng$se$model, nulls$se$null, tp$K, "se",
                                  macro, ztest)
        }
        res$R2_SE <- x$r2
        if (macro) res$R2_SE_macro <- x$r2_macro
        if (ztest && !is.null(x$ztest)) {
          res$Z_SE <- x$ztest$z
          res$pval_SE <- x$ztest$pval
        }
      }
      if ("chisq" %in% metrics) {
        x <- .optop_index_result_doc(eng$chisq$model, nulls$chisq$null,
                                     tp$K, "chisq", macro, ztest)
        res$R2_chisq <- x$r2
        if (macro) res$R2_chisq_macro <- x$r2_macro
        if (ztest && !is.null(x$ztest)) {
          res$Z_chisq <- x$ztest$z
          res$pval_chisq <- x$ztest$pval
        }
      }
      if ("deviance" %in% metrics) {
        x <- .optop_index_result_doc(eng$deviance$model, nulls$deviance$null,
                                     tp$K, "deviance", macro, ztest)
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
    if (is.null(block_size)) {
      block_size <- max(100L, min(W, floor(5e8 / (J * 8))))
    }
    acc <- matrix(0, W, 6)
    for (start in seq(1L, W, by = block_size)) {
      end <- min(start + block_size - 1L, W)
      w_idx <- start:end
      E_block <- if (do_model) {
        (theta %*% phi[, w_idx, drop = FALSE]) * partition$L
      } else {
        matrix(0, 0, 0)
      }
      acc[w_idx, ] <- optop_index_word_core(E_block, N, start - 1L,
                                            length(w_idx), L,
                                            pi_num[w_idx], eps,
                                            do_model, do_null,
                                            do_se, do_chisq, do_dev,
                                            n_threads)
    }
  } else {
    block_size_doc <- max(100L, min(J, floor(5e8 / (W * 8))))
    N_t <- Matrix::t(N)
    mask_bits <- optop_pack_mask_core(partition$rare_mask)
    tww <- if (do_model) t(phi) else matrix(0, 0, 0)
    acc <- matrix(0, J, 6)
    for (start in seq(1L, J, by = block_size_doc)) {
      end <- min(start + block_size_doc - 1L, J)
      j_idx <- start:end
      theta_blk <- if (do_model) {
        theta[j_idx, , drop = FALSE]
      } else {
        matrix(0, 0, 0)
      }
      acc[j_idx, ] <- optop_index_doc_core(tww, theta_blk, N_t, start - 1L,
                                           mask_bits, L[j_idx], pi_num, eps,
                                           do_model, do_null,
                                           do_se, do_chisq, do_dev,
                                           n_threads)
    }
  }
  list(se = list(model = acc[, 1], null = acc[, 2]),
       chisq = list(model = acc[, 3], null = acc[, 4]),
       deviance = list(model = acc[, 5], null = acc[, 6]))
}

# Assemble the document-level result list from the per-document
# discrepancies, exactly as the pre-engine implementations did.
.optop_index_result_doc <- function(D_K, D_null, K, metric, macro, ztest) {
  r2_doc <- ifelse(D_null > 0, 1 - D_K / D_null, 0)
  r2_micro <- 1 - sum(D_K) / sum(D_null)
  r2_macro <- mean(r2_doc[D_null > 0])
  result <- list(r2 = r2_micro, r2_macro = if (macro) r2_macro else NULL,
                 r2_doc = r2_doc, K = K, metric = metric)
  if (ztest) {
    result$ztest <- .optop_ztest(r2_doc[D_null > 0], r2_macro)
  }
  result
}

# Assemble the word-level result list from the per-word discrepancies.
.optop_index_result_word <- function(d_w, d_null, vocab, K, metric) {
  r2_word <- ifelse(d_null > 0, 1 - d_w / d_null, 0)
  names(r2_word) <- vocab
  list(r2_word = r2_word,
       r2_micro_word = sum((d_null / sum(d_null)) * r2_word),
       r2_macro_word = mean(r2_word[d_null > 0]),
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
  if (!identical(colnames(dtm), vocab_model))
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

  # partition alignment
  if (is.null(colnames(partition$rare_mask)) ||
      !identical(colnames(partition$rare_mask), vocab_model)) {
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
