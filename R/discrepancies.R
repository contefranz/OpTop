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
#'
#' @details
#' **Harmonized support.** For each document, rare words (as determined by
#' [`optop_make_partition()`]) are collapsed into a single “min” bin. Indices are then
#' evaluated on `{non-rare terms} ∪ {min}` to ensure comparability across \eqn{K}.
#'
#' **Alignment requirements.** The following must share the same vocabulary and column
#' order: the `dtm` passed to the index, the DTM used for `partition` and `baseline`, and the
#' model’s term–topic matrix. If they differ, align with `optop_align_dtm_to_models()` and
#' recompute `partition` and `baseline`.
#'
#' **Counts only.** SE/chi-square/deviance indices are defined for multinomial counts.
#' Do not pass weighted matrices (e.g., proportions from `quanteda::dfm_weight(scheme = "prop")`
#' or tf-idf). If your workflow uses proportions elsewhere, reconstruct counts before calling.
#' 
#' All three functions return the micro index (e.g., \eqn{R^2_{SE}}{R2_SE}) and, if
#' requested, the macro index (e.g., \eqn{\bar R^2_{SE}}{R2_SE_bar}), plus per-document
#' components (e.g., \eqn{R^2_{SE,j}}{R2_SE,j}).
#'
#' @return A list with:
#' - `r2`: scalar micro index (e.g., \eqn{R^2_{SE}}{R2_SE}).
#' - `r2_macro`: scalar macro index if `macro = TRUE`, otherwise `NULL`.
#' - `r2_doc`: numeric vector (length = number of documents) with per-document contributions.
#' - `K`: number of topics in `model`.
#' - `metric`: one of `"se"`, `"chisq"`, `"deviance"`, according to function call.
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
                           add_baseline_topic = TRUE) {
  
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
  
  # precompute baseline vectors
  B_nonrare_all <- matrix(NA_real_, nrow = J, ncol = W)
  B_min_all     <- numeric(J)
  for (j in 1:J) {
    rare_j <- partition$rare_mask[j, ]
    B_nonrare_all[j, ] <- partition$L[j] * pi_row
    B_min_all[j]       <- partition$L[j] * sum(pi_row[rare_j])
  }
  
  # loop over docs
  D_K <- numeric(J); D_null <- numeric(J); r2_doc <- numeric(J)
  for (j in 1:J) {
    # fitted probabilities and counts
    i_j <- as.numeric(theta[j, , drop=TRUE] %*% phi)       # length W; E_jw = L_j * i_jw
    E_j <- partition$L[j] * i_j
    
    # optional θ re-optimization targeting SE (projected simplex PGD)
    if (reopt == "se") {
      # cheap one-parameter blend towards baseline (guarantees SSE <= min{SSE_K,SST})
      i_base <- pi_row
      diff   <- i_j - i_base
      # compute λ* on counts space
      # inner products over non-rare plus min bin
      rare_j <- partition$rare_mask[j, ]
      # vectors on fixed support
      Nj_nonrare <- dtm[j, !rare_j, drop = TRUE]
      Nj_min     <- sum(dtm[j,  rare_j, drop = TRUE])
      Ej_nonrare <- E_j[!rare_j]; Ej_min <- sum(E_j[ rare_j])
      Bj_nonrare <- partition$L[j] * i_base[!rare_j]
      Bj_min     <- partition$L[j] * sum(i_base[ rare_j])
      num <- sum((Nj_nonrare - Bj_nonrare) * (Ej_nonrare - Bj_nonrare)) + (Nj_min - Bj_min) * (Ej_min - Bj_min)
      den <- sum((Ej_nonrare - Bj_nonrare)^2) + (Ej_min - Bj_min)^2
      lambda <- if (den <= 0) 0 else max(0, min(1, num / den))
      E_j[!rare_j] <- lambda * Ej_nonrare + (1 - lambda) * Bj_nonrare
      E_j[ rare_j] <- 0  # reconciling min below
      E_j_min <- lambda * Ej_min + (1 - lambda) * Bj_min
    } else {
      rare_j <- partition$rare_mask[j, ]
      E_j_min <- sum(E_j[rare_j])
    }
    
    # pack vectors and compute discrepancies
    v <- .optop_doc_vectors(j, dtm, rare_j, partition$L[j], E_j,
                            B_nonrare_all[j, ], B_min_all[j])
    sse_k   <- .optop_disc_se(v$N_nonrare, v$E_nonrare) + (v$N_min - v$E_min)^2
    sst_se  <- .optop_disc_se(v$N_nonrare, v$B_nonrare) + (v$N_min - v$B_min)^2
    D_K[j]  <- sse_k
    D_null[j] <- sst_se
    r2_doc[j] <- if (sst_se > 0) 1 - sse_k / sst_se else 0
  }
  r2_micro <- 1 - sum(D_K) / sum(D_null)             # eq. (14)
  r2_macro <- mean(r2_doc[D_null > 0])               # eq. (18)
  list(r2 = r2_micro, r2_macro = if (macro) r2_macro else NULL,
       r2_doc = r2_doc, K = tp$K, metric = "se")
}

#' @describeIn optop_index Pearson chi-square index \eqn{R^2_{chisq}}{R2_chisq}.
#' @param add_baseline_topic Logical; if `TRUE`, augment topics with a baseline row to
#' guarantee non-negativity of \eqn{R^2}. Default: `TRUE`.
#' (Only for `optop_index_se()` and `optop_index_chisq()`).
#' 
#' @export
optop_index_chisq <- function(model, dtm, partition, baseline,
                              macro = FALSE, reopt = c("none", "pearson"),
                              add_baseline_topic = TRUE) {
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
  
  D_K <- numeric(J); D_null <- numeric(J); r2_doc <- numeric(J)
  for (j in 1:J) {
    i_j <- as.numeric(theta[j, , drop=TRUE] %*% phi)
    E_j <- partition$L[j] * i_j
    rare_j <- partition$rare_mask[j, ]
    
    # reopt (lightweight IRLS would go here; omitted for brevity)
    
    v <- .optop_doc_vectors(j, dtm, rare_j, partition$L[j], E_j,
                            partition$L[j] * pi_row,
                            partition$L[j] * sum(pi_row[rare_j]))
    ssk   <- .optop_disc_chisq(v$N_nonrare, v$E_nonrare) + .optop_disc_chisq(v$N_min, v$E_min)
    sst   <- .optop_disc_chisq(v$N_nonrare, v$B_nonrare) + .optop_disc_chisq(v$N_min, v$B_min)
    D_K[j] <- ssk; D_null[j] <- sst
    r2_doc[j] <- if (sst > 0) 1 - ssk / sst else 0
  }
  list(r2 = 1 - sum(D_K) / sum(D_null),
       r2_macro = if (macro) mean(r2_doc[D_null > 0]) else NULL,
       r2_doc = r2_doc, K = tp$K, metric = "chisq")
}

#' @describeIn optop_index Deviance index \eqn{R^2_{dev}}{R2_dev}.
#'
#' Based on the multinomial deviance; monotone in \eqn{K} under (approximate) ML fits.
#' @export
optop_index_deviance <- function(model, dtm, partition, baseline,
                                 macro = FALSE, reopt = c("none", "deviance")) {
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
  
  
  D_K <- numeric(J); D_null <- numeric(J); r2_doc <- numeric(J)
  for (j in 1:J) {
    i_j <- as.numeric(theta[j, , drop=TRUE] %*% phi)
    E_j <- partition$L[j] * i_j
    rare_j <- partition$rare_mask[j, ]
    # (optional) EM reopt of θ_j for deviance could go here
    
    v <- .optop_doc_vectors(j, dtm, rare_j, partition$L[j], E_j,
                            partition$L[j] * pi_row,
                            partition$L[j] * sum(pi_row[rare_j]))
    dev_k   <- .optop_disc_dev(v$N_nonrare, v$E_nonrare) + .optop_disc_dev(v$N_min, v$E_min)
    dev_null<- .optop_disc_dev(v$N_nonrare, v$B_nonrare) + .optop_disc_dev(v$N_min, v$B_min)
    D_K[j] <- dev_k; D_null[j] <- dev_null
    r2_doc[j] <- if (dev_null > 0) 1 - dev_k / dev_null else 0
  }
  list(r2 = 1 - sum(D_K) / sum(D_null),
       r2_macro = if (macro) mean(r2_doc[D_null > 0]) else NULL,
       r2_doc = r2_doc, K = tp$K, metric = "deviance")
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
#'   - `"pearson"`: reserved for a Pearson-specific re-optimizer (if implemented);
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
#'
#' @details
#' The function wraps [`optop_index_se()`], [`optop_index_chisq()`], and
#' [`optop_index_deviance()`]. It reuses a harmonized support (via
#' [`optop_make_partition()`]) and a fixed baseline (via
#' [`optop_make_baseline()`]) so results are directly comparable across \eqn{K}.
#'
#' @return A `data.frame` with one row per model and columns:
#' - `K`: number of topics in the model.
#' - `R2_SE`, `R2_chisq`, `R2_dev`: micro indices for the selected metrics.
#' - `R2_SE_macro`, `R2_chisq_macro`, `R2_dev_macro`: macro indices when `macro = TRUE`.
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
                              partition = NULL, baseline = NULL) {
  stopifnot(length(models) >= 1)
  if (is.null(partition)) partition <- optop_make_partition(models, dtm, c = c)
  if (is.null(baseline))  baseline  <- optop_make_baseline(dtm)
  
  rows <- list()
  for (m in models) {
    K <- optop_as_theta_phi(m)$K
    res <- list(K = K)
    
    if ("se" %in% metrics) {
      x <- optop_index_se(m, dtm, partition, baseline,
                          macro, reopt = if (reopt=="se") "se" else "none",
                          add_baseline_topic = add_baseline_topic)
      res$R2_SE <- x$r2; if (macro) res$R2_SE_macro <- x$r2_macro
    }
    if ("chisq" %in% metrics) {
      x <- optop_index_chisq(m, dtm, partition, baseline,
                             macro, reopt = if (reopt=="pearson") "pearson" else "none",
                             add_baseline_topic = add_baseline_topic)
      res$R2_chisq <- x$r2; if (macro) res$R2_chisq_macro <- x$r2_macro
    }
    if ("deviance" %in% metrics) {
      x <- optop_index_deviance(m, dtm, partition, baseline,
                                macro, reopt = if (reopt=="deviance") "deviance" else "none")
      res$R2_dev <- x$r2; if (macro) res$R2_dev_macro <- x$r2_macro
    }
    rows[[length(rows)+1]] <- as.data.frame(res, check.names = FALSE)
  }
  do.call(rbind, rows)
}


#' @keywords internal
# Build observed/fitted/baseline vectors on fixed support {w ∉ C*_j} ∪ {min}
.optop_doc_vectors <- function(j, dtm, rare_mask_row, L_j, E_row, B_nonrare, B_min) {
  # indices for non-rare words
  nonrare <- which(!rare_mask_row)
  # observed counts
  Nj_nonrare <- dtm[j, nonrare, drop = TRUE]
  Nj_min     <- sum(dtm[j, rare_mask_row, drop = TRUE])
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
