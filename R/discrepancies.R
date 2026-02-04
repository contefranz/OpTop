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
#' @name optop_index
#' @aliases optop_index_se optop_index_chisq optop_index_deviance
NULL
#' @describeIn optop_index Squared-error index \eqn{R^2_{SE}}{R2_SE}.
#' @param add_baseline_topic Logical; if `TRUE`, augment topics with a baseline row to
#'   guarantee non-negativity of \eqn{R^2}. Default: `TRUE`.
#'   (Only for `optop_index_se()` and `optop_index_chisq()`).
#' @export
optop_index_se <- function(model, dtm, partition, baseline,
                           macro = FALSE, reopt = c("none", "se"),
                           add_baseline_topic = TRUE,
                           level = c("document", "word"),
                           block_size = NULL,
                           ztest = FALSE) {

  level <- match.arg(level)
  tp <- optop_as_theta_phi(model)
  vocab_model <- colnames(tp$phi)

  if (!identical(colnames(dtm), vocab_model))
    stop("DTM vocabulary/order differs from model. Use optop_align_dtm_to_models() and recompute partition/baseline.")
  # baseline alignment
  if (!is.null(names(baseline$pi_glob))) {
    pi_row <- baseline$pi_glob[vocab_model]
    if (any(is.na(pi_row)))
      stop("Baseline vocabulary does not match model. Recompute optop_make_baseline() on an aligned DTM.")
  } else {
    if (length(baseline$pi_glob) != ncol(tp$phi))
      stop("Baseline length != model vocab. Recompute baseline on aligned DTM.")
    pi_row <- baseline$pi_glob
  }

  # partition alignment
  if (is.null(colnames(partition$rare_mask)) ||
      !identical(colnames(partition$rare_mask), vocab_model)) {
    stop("Partition vocabulary != model vocabulary. Recompute optop_make_partition() on the aligned DTM.")
  }

  reopt <- match.arg(reopt)
  theta <- tp$theta
  phi <- tp$phi
  J <- nrow(theta)
  W <- ncol(phi)

  # optional baseline topic augmentation (guarantees R2 >= 0 for SE)
  if (add_baseline_topic) {
    phi <- rbind(phi, pi_row)
    # keep θ as-is; any θ-reopt may allocate weight to baseline topic
    theta <- cbind(theta, rep(0, J))
  }

  # =========================================================================
  # WORD-LEVEL AGGREGATION (block-based for memory efficiency)
  # =========================================================================
  if (level == "word") {
    # Auto-detect block size if not specified
    # Target: ~500 MB per block matrix (1.5 GB total for 3 matrices)
    if (is.null(block_size)) {
      block_size <- max(100L, min(W, floor(5e8 / (J * 8))))
    }

    SSE_w <- numeric(W)
    SST_w <- numeric(W)

    for (start in seq(1L, W, by = block_size)) {
      end <- min(start + block_size - 1L, W)
      w_idx <- start:end

      # Extract block of observed counts from sparse DTM
      # Converting only the block to dense: J × block_size
      N_block <- as.matrix(dtm[, w_idx, drop = FALSE])

      # Compute expected counts for this block: E_jw = L_j * Σ_k θ_jk φ_kw
      E_block <- (theta %*% phi[, w_idx, drop = FALSE]) * partition$L

      # Baseline expected counts: B_jw = L_j * π_glob(w)
      B_block <- outer(partition$L, pi_row[w_idx])

      # Squared-error discrepancies for this block
      SSE_w[w_idx] <- colSums((N_block - E_block)^2)
      SST_w[w_idx] <- colSums((N_block - B_block)^2)
    }

    # R² for each word
    r2_word <- ifelse(SST_w > 0, 1 - SSE_w / SST_w, 0)
    names(r2_word) <- vocab_model

    # Micro-Word: frequency-weighted average
    omega_w <- SST_w / sum(SST_w)
    r2_micro_word <- sum(omega_w * r2_word)

    # Macro-Word: unweighted average (excluding degenerate words)
    valid_words <- SST_w > 0
    r2_macro_word <- mean(r2_word[valid_words])

    return(list(r2_word = r2_word, r2_micro_word = r2_micro_word,
                r2_macro_word = r2_macro_word, K = tp$K, metric = "se"))
  }

  # =========================================================================
  # DOCUMENT-LEVEL AGGREGATION (block-based vectorization)
  # =========================================================================

  # Adaptive block size for documents: target ~500 MB per block matrix
  block_size_doc <- max(100L, min(J, floor(5e8 / (W * 8))))

  # Initialize output vectors

  D_K <- numeric(J)
  D_null <- numeric(J)
  r2_doc <- numeric(J)

  # Precompute B_min for all documents (vectorized)
  # B_min[j] = L[j] * sum(pi_row[rare_j])
  rare_mask <- partition$rare_mask  # J × W logical matrix
  B_min_all <- partition$L * as.numeric(rare_mask %*% pi_row)

  for (start in seq(1L, J, by = block_size_doc)) {
    end <- min(start + block_size_doc - 1L, J)
    j_idx <- start:end
    block_len <- length(j_idx)

    # Extract blocks: convert sparse DTM to dense for this block only
    N_block <- as.matrix(dtm[j_idx, , drop = FALSE])           # block × W
    L_block <- partition$L[j_idx]                              # block vector
    rare_block <- rare_mask[j_idx, , drop = FALSE]             # block × W

    # Compute expected counts: E_jw = L_j * Σ_k θ_jk φ_kw
    E_block <- (theta[j_idx, , drop = FALSE] %*% phi) * L_block  # block × W

    # Baseline expected counts: B_jw = L_j * π_glob(w)
    B_block <- outer(L_block, pi_row)                          # block × W

    if (reopt == "se") {
      # Re-optimization requires per-document lambda computation
      # Keep as nested loop within the block
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
    } else {
      # Vectorized computation using decomposition trick:
      # D_K = SSE_full - SSE_rare + SSE_min
      # where SSE_full = Σ_w (N-E)², SSE_rare = Σ_{rare w} (N-E)², SSE_min = (N_min - E_min)²

      diff_E <- N_block - E_block  # block × W
      diff_B <- N_block - B_block  # block × W

      # Model discrepancy D_K
      SSE_full <- rowSums(diff_E^2)
      SSE_rare <- rowSums(rare_block * diff_E^2)
      N_min <- rowSums(rare_block * N_block)
      E_min <- rowSums(rare_block * E_block)
      D_K[j_idx] <- SSE_full - SSE_rare + (N_min - E_min)^2

      # Baseline discrepancy D_null
      SST_full <- rowSums(diff_B^2)
      SST_rare <- rowSums(rare_block * diff_B^2)
      B_min <- B_min_all[j_idx]
      D_null[j_idx] <- SST_full - SST_rare + (N_min - B_min)^2

      # Per-document R²
      r2_doc[j_idx] <- ifelse(D_null[j_idx] > 0, 1 - D_K[j_idx] / D_null[j_idx], 0)
    }
  }

  r2_micro <- 1 - sum(D_K) / sum(D_null)
  r2_macro <- mean(r2_doc[D_null > 0])

  # Build result
  result <- list(r2 = r2_micro, r2_macro = if (macro) r2_macro else NULL,
                 r2_doc = r2_doc, K = tp$K, metric = "se")

  # Z-test for cross-document inference (Section 3.6.1)
  if (ztest) {
    result$ztest <- .optop_ztest(r2_doc[D_null > 0], r2_macro)
  }

  result
}

#' @describeIn optop_index Pearson chi-square index \eqn{R^2_{chisq}}{R2_chisq}.
#' @param add_baseline_topic Logical; if `TRUE`, augment topics with a baseline row to
#' guarantee non-negativity of \eqn{R^2}. Default: `TRUE`.
#' (Only for `optop_index_se()` and `optop_index_chisq()`).
#'
#' @export
optop_index_chisq <- function(model, dtm, partition, baseline,
                              macro = FALSE, reopt = c("none", "pearson"),
                              add_baseline_topic = TRUE,
                              level = c("document", "word"),
                              block_size = NULL,
                              ztest = FALSE) {

  level <- match.arg(level)
  tp <- optop_as_theta_phi(model)
  vocab_model <- colnames(tp$phi)

  if (!identical(colnames(dtm), vocab_model))
    stop("DTM vocabulary/order differs from model. Use optop_align_dtm_to_models() and recompute partition/baseline.")
  # baseline alignment
  if (!is.null(names(baseline$pi_glob))) {
    pi_row <- baseline$pi_glob[vocab_model]
    if (any(is.na(pi_row)))
      stop("Baseline vocabulary does not match model. Recompute optop_make_baseline() on an aligned DTM.")
  } else {
    if (length(baseline$pi_glob) != ncol(tp$phi))
      stop("Baseline length != model vocab. Recompute baseline on aligned DTM.")
    pi_row <- baseline$pi_glob
  }

  # partition alignment
  if (is.null(colnames(partition$rare_mask)) ||
      !identical(colnames(partition$rare_mask), vocab_model)) {
    stop("Partition vocabulary != model vocabulary. Recompute optop_make_partition() on the aligned DTM.")
  }
  reopt <- match.arg(reopt)
  theta <- tp$theta
  phi <- tp$phi
  J <- nrow(theta); W <- ncol(phi)

  if (add_baseline_topic) {
    phi <- rbind(phi, pi_row)
    theta <- cbind(theta, rep(0, J))
  }

  # =========================================================================
  # WORD-LEVEL AGGREGATION (block-based for memory efficiency)
  # =========================================================================
  if (level == "word") {
    eps <- 1e-12

    # Auto-detect block size if not specified
    if (is.null(block_size)) {
      block_size <- max(100L, min(W, floor(5e8 / (J * 8))))
    }

    chisq_w <- numeric(W)
    chisq_w_null <- numeric(W)

    for (start in seq(1L, W, by = block_size)) {
      end <- min(start + block_size - 1L, W)
      w_idx <- start:end

      # Extract block of observed counts from sparse DTM
      N_block <- as.matrix(dtm[, w_idx, drop = FALSE])

      # Compute expected counts for this block
      E_block <- (theta %*% phi[, w_idx, drop = FALSE]) * partition$L
      E_block <- pmax(E_block, eps)

      # Baseline expected counts
      B_block <- outer(partition$L, pi_row[w_idx])
      B_block <- pmax(B_block, eps)

      # Chi-square discrepancies for this block
      chisq_w[w_idx] <- colSums((N_block - E_block)^2 / E_block)
      chisq_w_null[w_idx] <- colSums((N_block - B_block)^2 / B_block)
    }

    # R² for each word
    r2_word <- ifelse(chisq_w_null > 0, 1 - chisq_w / chisq_w_null, 0)
    names(r2_word) <- vocab_model

    # Micro-Word: frequency-weighted average
    omega_w <- chisq_w_null / sum(chisq_w_null)
    r2_micro_word <- sum(omega_w * r2_word)

    # Macro-Word: unweighted average (excluding degenerate words)
    valid_words <- chisq_w_null > 0
    r2_macro_word <- mean(r2_word[valid_words])

    return(list(r2_word = r2_word, r2_micro_word = r2_micro_word,
                r2_macro_word = r2_macro_word, K = tp$K, metric = "chisq"))
  }

  # =========================================================================
  # DOCUMENT-LEVEL AGGREGATION (block-based vectorization)
  # =========================================================================

  eps <- 1e-12

  # Adaptive block size for documents: target ~500 MB per block matrix
  block_size_doc <- max(100L, min(J, floor(5e8 / (W * 8))))

  # Initialize output vectors
  D_K <- numeric(J)
  D_null <- numeric(J)
  r2_doc <- numeric(J)

  # Precompute B_min for all documents (vectorized)
  rare_mask <- partition$rare_mask  # J × W logical matrix
  B_min_all <- partition$L * as.numeric(rare_mask %*% pi_row)

  for (start in seq(1L, J, by = block_size_doc)) {
    end <- min(start + block_size_doc - 1L, J)
    j_idx <- start:end

    # Extract blocks: convert sparse DTM to dense for this block only
    N_block <- as.matrix(dtm[j_idx, , drop = FALSE])           # block × W
    L_block <- partition$L[j_idx]                              # block vector
    rare_block <- rare_mask[j_idx, , drop = FALSE]             # block × W

    # Compute expected counts: E_jw = L_j * Σ_k θ_jk φ_kw
    E_block <- (theta[j_idx, , drop = FALSE] %*% phi) * L_block  # block × W
    E_block <- pmax(E_block, eps)

    # Baseline expected counts: B_jw = L_j * π_glob(w)
    B_block <- outer(L_block, pi_row)                          # block × W
    B_block <- pmax(B_block, eps)

    # Vectorized chi-square computation using decomposition:
    # D_K = χ²_full - χ²_rare + χ²_min

    # Model discrepancy D_K
    chisq_full <- rowSums((N_block - E_block)^2 / E_block)
    chisq_rare <- rowSums(rare_block * (N_block - E_block)^2 / E_block)
    N_min <- rowSums(rare_block * N_block)
    E_min <- pmax(rowSums(rare_block * E_block), eps)
    chisq_min <- (N_min - E_min)^2 / E_min
    D_K[j_idx] <- chisq_full - chisq_rare + chisq_min

    # Baseline discrepancy D_null
    chisq_full_null <- rowSums((N_block - B_block)^2 / B_block)
    chisq_rare_null <- rowSums(rare_block * (N_block - B_block)^2 / B_block)
    B_min <- pmax(B_min_all[j_idx], eps)
    chisq_min_null <- (N_min - B_min)^2 / B_min
    D_null[j_idx] <- chisq_full_null - chisq_rare_null + chisq_min_null

    # Per-document R²
    r2_doc[j_idx] <- ifelse(D_null[j_idx] > 0, 1 - D_K[j_idx] / D_null[j_idx], 0)
  }

  r2_micro <- 1 - sum(D_K) / sum(D_null)
  r2_macro <- mean(r2_doc[D_null > 0])

  # Build result
  result <- list(r2 = r2_micro,
                 r2_macro = if (macro) r2_macro else NULL,
                 r2_doc = r2_doc, K = tp$K, metric = "chisq")

  # Z-test for cross-document inference (Section 3.6.1)
  if (ztest) {
    result$ztest <- .optop_ztest(r2_doc[D_null > 0], r2_macro)
  }

  result
}

#' @describeIn optop_index Deviance index \eqn{R^2_{dev}}{R2_dev}.
#'
#' Based on the multinomial deviance; monotone in \eqn{K} under (approximate) ML fits.
#' @export
optop_index_deviance <- function(model, dtm, partition, baseline,
                                 macro = FALSE, reopt = c("none", "deviance"),
                                 level = c("document", "word"),
                                 block_size = NULL,
                                 ztest = FALSE) {

  level <- match.arg(level)
  tp <- optop_as_theta_phi(model)
  vocab_model <- colnames(tp$phi)

  if (!identical(colnames(dtm), vocab_model))
    stop("DTM vocabulary/order differs from model. Use optop_align_dtm_to_models() and recompute partition/baseline.")
  # baseline alignment
  if (!is.null(names(baseline$pi_glob))) {
    pi_row <- baseline$pi_glob[vocab_model]
    if (any(is.na(pi_row)))
      stop("Baseline vocabulary does not match model. Recompute optop_make_baseline() on an aligned DTM.")
  } else {
    if (length(baseline$pi_glob) != ncol(tp$phi))
      stop("Baseline length != model vocab. Recompute baseline on aligned DTM.")
    pi_row <- baseline$pi_glob
  }

  # partition alignment
  if (is.null(colnames(partition$rare_mask)) ||
      !identical(colnames(partition$rare_mask), vocab_model)) {
    stop("Partition vocabulary != model vocabulary. Recompute optop_make_partition() on the aligned DTM.")
  }
  reopt <- match.arg(reopt)
  theta <- tp$theta
  phi <- tp$phi
  J <- nrow(theta); W <- ncol(phi)

  # =========================================================================
  # WORD-LEVEL AGGREGATION (block-based for memory efficiency)
  # =========================================================================
  if (level == "word") {
    eps <- 1e-12

    # Auto-detect block size if not specified
    if (is.null(block_size)) {
      block_size <- max(100L, min(W, floor(5e8 / (J * 8))))
    }

    dev_w <- numeric(W)
    dev_w_null <- numeric(W)

    for (start in seq(1L, W, by = block_size)) {
      end <- min(start + block_size - 1L, W)
      w_idx <- start:end

      # Extract block of observed counts from sparse DTM
      N_block <- as.matrix(dtm[, w_idx, drop = FALSE])

      # Compute expected counts for this block
      E_block <- (theta %*% phi[, w_idx, drop = FALSE]) * partition$L
      E_block <- pmax(E_block, eps)

      # Baseline expected counts
      B_block <- outer(partition$L, pi_row[w_idx])
      B_block <- pmax(B_block, eps)

      # Deviance for each word in block (only where N > 0)
      for (i in seq_along(w_idx)) {
        w <- w_idx[i]
        N_w <- N_block[, i]
        E_w <- E_block[, i]
        B_w <- B_block[, i]
        idx <- N_w > 0
        if (any(idx)) {
          dev_w[w] <- 2 * sum(N_w[idx] * (log(N_w[idx]) - log(E_w[idx])))
          dev_w_null[w] <- 2 * sum(N_w[idx] * (log(N_w[idx]) - log(B_w[idx])))
        }
      }
    }

    # R² for each word
    r2_word <- ifelse(dev_w_null > 0, 1 - dev_w / dev_w_null, 0)
    names(r2_word) <- vocab_model

    # Micro-Word: frequency-weighted average
    omega_w <- dev_w_null / sum(dev_w_null)
    r2_micro_word <- sum(omega_w * r2_word)

    # Macro-Word: unweighted average (excluding degenerate words)
    valid_words <- dev_w_null > 0
    r2_macro_word <- mean(r2_word[valid_words])

    return(list(r2_word = r2_word, r2_micro_word = r2_micro_word,
                r2_macro_word = r2_macro_word, K = tp$K, metric = "deviance"))
  }

  # =========================================================================
  # DOCUMENT-LEVEL AGGREGATION (block-based vectorization)
  # =========================================================================

  eps <- 1e-12

  # Adaptive block size for documents: target ~500 MB per block matrix
  block_size_doc <- max(100L, min(J, floor(5e8 / (W * 8))))

  # Initialize output vectors
  D_K <- numeric(J)
  D_null <- numeric(J)
  r2_doc <- numeric(J)

  # Precompute B_min for all documents (vectorized)
  rare_mask <- partition$rare_mask  # J × W logical matrix
  B_min_all <- partition$L * as.numeric(rare_mask %*% pi_row)

  # Helper: safe deviance contribution (handles 0*log(0) = 0 convention)
  safe_dev_contrib <- function(N, E) {
    E <- pmax(E, eps)
    result <- 2 * N * log(N / E)
    result[N == 0] <- 0
    result
  }

  for (start in seq(1L, J, by = block_size_doc)) {
    end <- min(start + block_size_doc - 1L, J)
    j_idx <- start:end

    # Extract blocks: convert sparse DTM to dense for this block only
    N_block <- as.matrix(dtm[j_idx, , drop = FALSE])           # block × W
    L_block <- partition$L[j_idx]                              # block vector
    rare_block <- rare_mask[j_idx, , drop = FALSE]             # block × W

    # Compute expected counts: E_jw = L_j * Σ_k θ_jk φ_kw
    E_block <- (theta[j_idx, , drop = FALSE] %*% phi) * L_block  # block × W

    # Baseline expected counts: B_jw = L_j * π_glob(w)
    B_block <- outer(L_block, pi_row)                          # block × W

    # Vectorized deviance computation using decomposition:
    # D_K = Dev_full - Dev_rare + Dev_min

    # Model discrepancy D_K
    dev_contrib_E <- safe_dev_contrib(N_block, E_block)        # block × W
    dev_full <- rowSums(dev_contrib_E)
    dev_rare <- rowSums(rare_block * dev_contrib_E)
    N_min <- rowSums(rare_block * N_block)
    E_min <- pmax(rowSums(rare_block * E_block), eps)
    dev_min <- 2 * N_min * log(pmax(N_min, eps) / E_min)
    dev_min[N_min == 0] <- 0
    D_K[j_idx] <- dev_full - dev_rare + dev_min

    # Baseline discrepancy D_null
    dev_contrib_B <- safe_dev_contrib(N_block, B_block)        # block × W
    dev_full_null <- rowSums(dev_contrib_B)
    dev_rare_null <- rowSums(rare_block * dev_contrib_B)
    B_min <- pmax(B_min_all[j_idx], eps)
    dev_min_null <- 2 * N_min * log(pmax(N_min, eps) / B_min)
    dev_min_null[N_min == 0] <- 0
    D_null[j_idx] <- dev_full_null - dev_rare_null + dev_min_null

    # Per-document R²
    r2_doc[j_idx] <- ifelse(D_null[j_idx] > 0, 1 - D_K[j_idx] / D_null[j_idx], 0)
  }

  r2_micro <- 1 - sum(D_K) / sum(D_null)
  r2_macro <- mean(r2_doc[D_null > 0])

  # Build result
  result <- list(r2 = r2_micro,
                 r2_macro = if (macro) r2_macro else NULL,
                 r2_doc = r2_doc, K = tp$K, metric = "deviance")

  # Z-test for cross-document inference (Section 3.6.1)
  if (ztest) {
    result$ztest <- .optop_ztest(r2_doc[D_null > 0], r2_macro)
  }

  result
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
#'   during word-level computation. If `NULL` (default), automatically chosen based on
#'   corpus size. Only used when `level = "word"`.
#' @param ztest Logical; if `TRUE`, append Z-test statistics for cross-document inference.
#'   Only applicable when `level = "document"`. Default: `FALSE`.
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
#' @export
optop_index_table <- function(models, dtm, metrics = c("se","chisq","deviance"),
                              c = 5, macro = FALSE, reopt = "none",
                              add_baseline_topic = TRUE,
                              partition = NULL, baseline = NULL,
                              level = c("document", "word"),
                              block_size = NULL,
                              ztest = FALSE) {

  level <- match.arg(level)
  stopifnot(length(models) >= 1)
  if (is.null(partition)) partition <- optop_make_partition(models, dtm, c = c)
  if (is.null(baseline))  baseline  <- optop_make_baseline(dtm)

  rows <- list()
  for (m in models) {
    K <- optop_as_theta_phi(m)$K
    res <- list(K = K)

    if (level == "word") {
      # Word-level indices
      if ("se" %in% metrics) {
        x <- optop_index_se(m, dtm, partition, baseline,
                            macro = FALSE, reopt = if (reopt=="se") "se" else "none",
                            add_baseline_topic = add_baseline_topic,
                            level = "word", block_size = block_size, ztest = FALSE)
        res$R2_SE_micro_word <- x$r2_micro_word
        res$R2_SE_macro_word <- x$r2_macro_word
      }
      if ("chisq" %in% metrics) {
        x <- optop_index_chisq(m, dtm, partition, baseline,
                               macro = FALSE, reopt = if (reopt=="pearson") "pearson" else "none",
                               add_baseline_topic = add_baseline_topic,
                               level = "word", block_size = block_size, ztest = FALSE)
        res$R2_chisq_micro_word <- x$r2_micro_word
        res$R2_chisq_macro_word <- x$r2_macro_word
      }
      if ("deviance" %in% metrics) {
        x <- optop_index_deviance(m, dtm, partition, baseline,
                                  macro = FALSE, reopt = if (reopt=="deviance") "deviance" else "none",
                                  level = "word", block_size = block_size, ztest = FALSE)
        res$R2_dev_micro_word <- x$r2_micro_word
        res$R2_dev_macro_word <- x$r2_macro_word
      }
    } else {
      # Document-level indices
      if ("se" %in% metrics) {
        x <- optop_index_se(m, dtm, partition, baseline,
                            macro, reopt = if (reopt=="se") "se" else "none",
                            add_baseline_topic = add_baseline_topic,
                            level = "document", block_size = block_size, ztest = ztest)
        res$R2_SE <- x$r2
        if (macro) res$R2_SE_macro <- x$r2_macro
        if (ztest && !is.null(x$ztest)) {
          res$Z_SE <- x$ztest$z
          res$pval_SE <- x$ztest$pval
        }
      }
      if ("chisq" %in% metrics) {
        x <- optop_index_chisq(m, dtm, partition, baseline,
                               macro, reopt = if (reopt=="pearson") "pearson" else "none",
                               add_baseline_topic = add_baseline_topic,
                               level = "document", block_size = block_size, ztest = ztest)
        res$R2_chisq <- x$r2
        if (macro) res$R2_chisq_macro <- x$r2_macro
        if (ztest && !is.null(x$ztest)) {
          res$Z_chisq <- x$ztest$z
          res$pval_chisq <- x$ztest$pval
        }
      }
      if ("deviance" %in% metrics) {
        x <- optop_index_deviance(m, dtm, partition, baseline,
                                  macro, reopt = if (reopt=="deviance") "deviance" else "none",
                                  level = "document", block_size = block_size, ztest = ztest)
        res$R2_dev <- x$r2
        if (macro) res$R2_dev_macro <- x$r2_macro
        if (ztest && !is.null(x$ztest)) {
          res$Z_dev <- x$ztest$z
          res$pval_dev <- x$ztest$pval
        }
      }
    }
    rows[[length(rows)+1]] <- as.data.frame(res, check.names = FALSE)
  }
  out <- do.call(rbind, rows)
  data.table::setDT(out)
}


#' @keywords internal
# Build observed/fitted/baseline vectors on fixed support {w ∉ C*_j} ∪ {min}
.optop_doc_vectors <- function(j, dtm, rare_mask_row, L_j, E_row, B_nonrare, B_min) {
  # indices for non-rare words
  nonrare <- which(!rare_mask_row)
  # observed counts (as.numeric ensures plain vector even from S4 dfm/dgCMatrix)
  Nj_nonrare <- as.numeric(dtm[j, nonrare, drop = TRUE])
  Nj_min     <- sum(as.numeric(dtm[j, rare_mask_row, drop = TRUE]))
  # fitted expected counts (already E_row = L_j * i_j)
  EK_nonrare <- E_row[nonrare]
  EK_min     <- sum(E_row[rare_mask_row])
  list(N_nonrare = Nj_nonrare, N_min = Nj_min,
       E_nonrare = EK_nonrare, E_min = EK_min,
       B_nonrare = B_nonrare[nonrare], B_min = B_min)
}

.optop_disc_se <- function(N, E) sum((N - E)^2)
.optop_disc_chisq <- function(N, E, eps = 1e-12) {
  E <- pmax(E, eps); sum((N - E)^2 / E)
}
.optop_disc_dev <- function(N, E, eps = 1e-12) {

  E <- pmax(E, eps); idx <- N > 0
  2 * sum(N[idx] * (log(N[idx]) - log(E[idx])))
}

#' Z-test for cross-document inference on Macro R² index
#'
#' Tests whether the K-topic model provides statistically significant
#' improvement over the no-topics baseline.
#'
#' @param r2_doc Numeric vector of document-level R² values (length J).
#' @param r2_macro Scalar macro index (mean of r2_doc over valid documents).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{z}: Z-statistic.
#'   \item \code{pval}: One-sided p-value for H1: μ_R² > 0.
#'   \item \code{se}: Standard error σ̂_R.
#'   \item \code{ci}: 95\% confidence interval for the true mean R².
#'   \item \code{J}: Number of valid documents used.
#' }
#'
#' @details
#' The null hypothesis is H0: μ_R² ≤ 0 (the topic model is no better than
#' the global distribution on average). Under regularity conditions, the
#' test statistic Z = √J · R̄²_Macro / σ̂_R is asymptotically N(0,1).
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
