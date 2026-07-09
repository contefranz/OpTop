#' Held-out goodness-of-fit indices with conditional inference
#'
#' Evaluate the OpTop goodness-of-fit indices on an independent evaluation
#' sample, following Section 3.7 of Lewis and Grossetti (2026). The topic
#' models and the global baseline are estimated on the training sample; for
#' each evaluation document the topic weights are folded in from that
#' document and the trained topics, the harmonized support is rebuilt on the
#' evaluation corpus with the training baseline as the null member, and the
#' discrepancies are evaluated on the evaluation documents. The document
#' scores are conditionally i.i.d. given the training sample, which yields
#' confidence intervals for the average held-out fit and paired comparisons
#' across topic counts ([optop_gain_table()]).
#'
#' @param models A list of topic models fitted on the \emph{training}
#'   corpus, spanning a grid of \eqn{K}; every class supported by the
#'   adapters and by fold-in (topicmodels fits, seededlda fits, NLPstudio
#'   `nlp_topic_fit` with a stored backend, text2vec WarpLDA via
#'   [optop_warplda()]).
#' @param dtm_eval A counts document-term matrix of the \emph{evaluation}
#'   documents. Columns are matched to the training vocabulary by name;
#'   evaluation-only words follow the out-of-support convention below.
#' @param baseline The \emph{training} baseline from
#'   [optop_make_baseline()], optionally smoothed (`smooth_lambda`).
#' @param c Positive scalar of the harmonized threshold
#'   \eqn{\tau_j = c / L_j}; default `1`.
#' @param metrics Character subset of `c("se", "chisq", "deviance")`.
#' @param conf Confidence level of the Macro interval (default `0.95`).
#' @param stabilize Nonnegative scalar \eqn{\varepsilon}{epsilon} of the
#'   stabilized index
#'   \eqn{1 - D^{ho}_j(K) / \max\{D^{ho}_j(\mathrm{null}), \varepsilon\}}{1 - D_j(K)/max(D_j(null), epsilon)}.
#'   Default `0` (no stabilization; degenerate documents carry `NA` and are
#'   excluded). A positive value makes every document score finite and
#'   changes the estimand to the mean of the stabilized score; fix it before
#'   evaluation and report sensitivity.
#' @param n_threads Integer; OpenMP threads of the compiled kernels.
#' @param ... Passed to the engine's fold-in routine (e.g. sampler
#'   iterations).
#'
#' @details
#' **Held-out target.** The scores implement held-out-document
#' reconstruction: the same evaluation document is used to infer its topic
#' mixture and to score its fit. For held-out-token completion, split each
#' evaluation document upstream and pass the fold-in tokens through the
#' engine directly.
#'
#' **Out-of-support convention.** The evaluation vocabulary is the training
#' support. Tokens of words absent from the training vocabulary are assigned
#' to the observed min-bin of their document (through a synthetic always-rare
#' column with zero fitted and baseline probability), and document lengths
#' \eqn{L_j} keep every token. Alternatively, smooth the training baseline
#' (`optop_make_baseline(dtm_train, smooth_lambda)`) and report sensitivity.
#'
#' **Inference.** For each \eqn{K}, the Macro index is the mean of the
#' document scores over the non-degenerate evaluation documents
#' \eqn{J_+}{J+}; its standard error uses the \eqn{1/(J_+ - 1)}{1/(J+ - 1)}
#' sample variance and the interval is
#' \eqn{\mathrm{Macro} \pm z_{1-\alpha/2}\, \hat\sigma / \sqrt{J_+}}{Macro +/- z * sigma_hat / sqrt(J+)}
#' (Proposition 2 of the paper). The Micro index and the Micro-Macro gap are
#' smooth functions of conditional means, and their standard errors follow
#' by the delta method. All inference is conditional on the training sample:
#' the target is the expected per-document fit under the training-fitted
#' global objects, not an unconditional population quantity.
#'
#' **Documents dropped.** Evaluation documents with zero aligned tokens
#' cannot be folded in and are dropped with a warning.
#'
#' @return An object of class `optop_holdout`: a list with
#' - `summary`: `data.table` with one row per (metric, K): `macro`,
#'   `macro_se`, `ci_lo`, `ci_hi`, `micro`, `micro_se`, `gap`, `gap_se`
#'   (gap = micro − macro), and `n_docs` (the size of \eqn{J_+}{J+}).
#' - `scores`: named list (one element per metric) of numeric matrices,
#'   evaluation documents by \eqn{K}, holding the per-document held-out
#'   indices (`NA` on degenerate documents when `stabilize = 0`).
#' - `d_model`, `d_null`: same shape, the per-document discrepancies.
#' - `K`, `metrics`, `conf`, `stabilize`, `c`: the evaluation design.
#'
#' @examples
#' \dontrun{
#' train <- dfm_counts[1:40, ]
#' eval <- dfm_counts[41:59, ]
#' models <- lapply(seq(2, 20, 2), function(k) topicmodels::LDA(train, k = k))
#' base_tr <- optop_make_baseline(train)
#' ho <- optop_index_holdout(models, eval, base_tr)
#' ho$summary
#' optop_gain_table(ho, epsilon = 0.01)
#' }
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2026). Goodness-of-fit indices and
#' diagnostics for topic models. Working paper.
#'
#' @seealso [optop_gain_table()], [optop_moment_test()],
#'   [optop_make_baseline()], [optop_index_table()]
#'
#' @export
optop_index_holdout <- function(models, dtm_eval, baseline, c = 1,
                                metrics = c("se", "chisq", "deviance"),
                                conf = 0.95, stabilize = 0,
                                n_threads = 1L, ...) {
  stopifnot(length(models) >= 1)
  metrics <- intersect(c("se", "chisq", "deviance"), metrics)
  if (!length(metrics)) stop("no valid metric requested")
  if (!is.numeric(conf) || length(conf) != 1L || conf <= 0 || conf >= 1) {
    stop("conf must be a single number in (0, 1)")
  }
  if (!is.numeric(stabilize) || length(stabilize) != 1L ||
      !is.finite(stabilize) || stabilize < 0) {
    stop("stabilize must be a single nonnegative number")
  }

  prep <- .optop_holdout_prepare(models, dtm_eval, baseline, ...)

  # held-out harmonized partition on the evaluation corpus: the union runs
  # over K0 = K U {null} with the TRAINING baseline as the null member
  part <- optop_make_partition(prep$tp_eval, prep$dtm_aug, c = c,
                               n_threads = n_threads,
                               pi_glob = prep$pi_aug)
  if ("chisq" %in% metrics) .optop_report_chisq_min(part)
  baseline_aug <- list(pi_glob = stats::setNames(prep$pi_aug,
                                                 prep$vocab_aug))
  pi_row <- .optop_validate_alignment(prep$vocab_aug, prep$dtm_aug, part,
                                      baseline_aug)

  # the baseline (no-topics) side is model-independent: one engine sweep
  nulls <- .optop_index_engine(theta = NULL, phi = NULL, prep$dtm_aug, part,
                               pi_row, metrics, "document", NULL, n_threads,
                               do_model = FALSE, do_null = TRUE)

  n_k <- length(models)
  J <- nrow(prep$dtm_aug)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  empty_mat <- function() {
    matrix(NA_real_, J, n_k,
           dimnames = list(rownames(prep$dtm_aug), prep$K))
  }
  scores <- d_model <- d_null <- stats::setNames(
    lapply(metrics, function(m) empty_mat()), metrics
  )
  rows <- list()

  for (i_mod in seq_len(n_k)) {
    tp <- prep$tp_eval[[i_mod]]
    eng <- .optop_index_engine(tp$theta, tp$phi, prep$dtm_aug, part, pi_row,
                               metrics, "document", NULL, n_threads,
                               do_model = TRUE, do_null = FALSE)
    for (metric in metrics) {
      D_K <- eng[[metric]]$model
      D_0 <- nulls[[metric]]$null
      st <- .optop_holdout_stats(D_K, D_0, stabilize, z)
      scores[[metric]][, i_mod] <- st$r
      d_model[[metric]][, i_mod] <- D_K
      d_null[[metric]][, i_mod] <- D_0
      rows[[length(rows) + 1L]] <- data.table::data.table(
        metric = metric, K = prep$K[i_mod],
        macro = st$macro, macro_se = st$macro_se,
        ci_lo = st$macro - z * st$macro_se,
        ci_hi = st$macro + z * st$macro_se,
        micro = st$micro, micro_se = st$micro_se,
        gap = st$micro - st$macro, gap_se = st$gap_se,
        n_docs = st$n
      )
    }
  }

  out <- list(summary = data.table::rbindlist(rows),
              scores = scores, d_model = d_model, d_null = d_null,
              K = prep$K, metrics = metrics, conf = conf,
              stabilize = stabilize, c = c)
  class(out) <- c("optop_holdout", "list")
  out
}

#' Shared held-out preparation
#'
#' Validates the training objects, aligns the evaluation counts to the
#' training vocabulary, applies the out-of-support convention, folds the
#' document-topic weights in for every model, and returns the augmented
#' evaluation objects the held-out tools consume.
#'
#' @param models,dtm_eval,baseline As in [optop_index_holdout()].
#' @param ... Passed to [optop_fold_in()].
#'
#' @return A list with the augmented dtm (`dtm_aug`), augmented baseline
#'   probabilities and vocabulary (`pi_aug`, `vocab_aug`), the per-model
#'   evaluation containers (`tp_eval`, theta folded in and phi augmented),
#'   the aligned unaugmented counts (`dtm_aligned`), the topic grid (`K`)
#'   and the training vocabulary (`vocab_tr`).
#'
#' @keywords internal
.optop_holdout_prepare <- function(models, dtm_eval, baseline, ...) {
  if (is.null(names(baseline$pi_glob))) {
    stop("baseline$pi_glob must be named; recompute optop_make_baseline()")
  }
  vocab_tr <- names(baseline$pi_glob)

  tp_tr <- lapply(models, optop_as_theta_phi)
  terms1 <- tp_tr[[1L]]$terms
  same <- vapply(tp_tr, function(tp) identical(tp$terms, terms1), logical(1))
  if (!all(same)) {
    stop("all models must share the same vocabulary in the same order")
  }
  if (!identical(terms1, vocab_tr)) {
    stop(paste("the training baseline vocabulary differs from the model",
               "vocabulary; recompute optop_make_baseline() on the aligned",
               "training dtm"))
  }
  K <- vapply(tp_tr, `[[`, numeric(1), "K")
  if (is.unsorted(K, strictly = TRUE)) {
    o <- order(K)
    if (anyDuplicated(K)) stop("one model per topic count is required")
    models <- models[o]
    tp_tr <- tp_tr[o]
    K <- K[o]
  }

  # align the evaluation counts to the training support; words the models
  # never saw follow the out-of-support convention (observed min-bin mass)
  dtm_eval <- methods::as(dtm_eval, "CsparseMatrix")
  if (is.null(colnames(dtm_eval))) {
    stop("dtm_eval must have column names to match the training vocabulary")
  }
  in_tr <- colnames(dtm_eval) %in% vocab_tr
  oov <- Matrix::rowSums(dtm_eval[, !in_tr, drop = FALSE])
  idx <- match(vocab_tr, colnames(dtm_eval))
  hit <- !is.na(idx)
  dtm_aligned <- Matrix::sparseMatrix(
    i = integer(0), j = integer(0), x = numeric(0),
    dims = c(nrow(dtm_eval), length(vocab_tr)),
    dimnames = list(rownames(dtm_eval), vocab_tr)
  )
  dtm_aligned[, which(hit)] <- dtm_eval[, idx[hit], drop = FALSE]

  # documents without any aligned token cannot be folded in
  keep <- Matrix::rowSums(dtm_aligned) > 0
  if (!all(keep)) {
    cli::cli_warn(paste(
      "dropped {sum(!keep)} evaluation document{?s} with no token on the",
      "training vocabulary"
    ))
    dtm_aligned <- dtm_aligned[keep, , drop = FALSE]
    oov <- oov[keep]
  }

  has_oov <- any(oov > 0)
  if (has_oov) {
    cli::cli_alert_info(paste(
      "{sum(oov > 0)} evaluation document{?s} carr{?ies/y} out-of-support",
      "words ({sprintf('%.4f', sum(oov) / (sum(dtm_aligned) + sum(oov)))}",
      "of the evaluation mass); their tokens enter the observed min-bin"
    ))
    # one synthetic column with zero fitted and baseline probability: it is
    # always rare (p = 0 < tau_j), so its counts land in N_min while E_min
    # and B_min are unchanged, and L_j keeps every token
    oov_col <- Matrix::Matrix(oov, ncol = 1, sparse = TRUE,
                              dimnames = list(rownames(dtm_aligned),
                                              "..optop_oov.."))
    dtm_aug <- methods::as(cbind(dtm_aligned, oov_col), "CsparseMatrix")
    vocab_aug <- c(vocab_tr, "..optop_oov..")
    pi_aug <- c(unname(baseline$pi_glob), 0)
  } else {
    dtm_aug <- dtm_aligned
    vocab_aug <- vocab_tr
    pi_aug <- unname(baseline$pi_glob)
  }

  # fold the document-topic weights in, one model at a time, and augment
  # phi with the zero out-of-support column where needed
  tp_eval <- vector("list", length(models))
  for (i_mod in seq_along(models)) {
    theta <- optop_fold_in(models[[i_mod]], dtm_aligned, ...)
    if (!is.matrix(theta) || nrow(theta) != nrow(dtm_aligned)) {
      stop("fold-in returned no valid document-topic matrix")
    }
    rownames(theta) <- rownames(dtm_aligned)
    phi <- tp_tr[[i_mod]]$phi
    if (has_oov) {
      phi <- cbind(phi, 0)
    }
    colnames(phi) <- vocab_aug
    tp_eval[[i_mod]] <- .optop_tp(theta, phi)
  }

  list(dtm_aug = dtm_aug, dtm_aligned = dtm_aligned,
       pi_aug = pi_aug, vocab_aug = vocab_aug, vocab_tr = vocab_tr,
       tp_eval = tp_eval, K = K)
}

#' Held-out score statistics for one (metric, K) cell
#'
#' Computes the per-document scores and the Macro, Micro, and gap estimates
#' with their standard errors: Macro from the i.i.d. sample variance over
#' J+, Micro = 1 - mean(D_K)/mean(D_null) and gap = micro - macro by the
#' delta method on the sample means of (D_K, D_null, r).
#'
#' @param D_K,D_null Per-document fitted and baseline discrepancies.
#' @param stabilize The stabilization constant (0 = none).
#' @param z The normal quantile of the confidence level (unused here beyond
#'   symmetry with the caller; kept for clarity).
#'
#' @return A list with `r`, `macro`, `macro_se`, `micro`, `micro_se`,
#'   `gap_se`, and `n` (the J+ count).
#'
#' @keywords internal
.optop_holdout_stats <- function(D_K, D_null, stabilize, z) {
  if (stabilize > 0) {
    # stabilized score: every document is finite, the estimand is the mean
    # of the stabilized score
    r <- 1 - D_K / pmax(D_null, stabilize)
    valid <- rep(TRUE, length(r))
  } else {
    valid <- D_null > 0
    r <- rep(NA_real_, length(D_K))
    r[valid] <- 1 - D_K[valid] / D_null[valid]
  }
  n <- sum(valid)
  rv <- r[valid]
  A <- D_K[valid]
  B <- D_null[valid]

  macro <- mean(rv)
  macro_se <- stats::sd(rv) / sqrt(n)

  # micro = 1 - Abar/Bbar; delta method with gradient (-1/Bbar, Abar/Bbar^2)
  Abar <- mean(A)
  Bbar <- mean(B)
  micro <- 1 - Abar / Bbar
  g2 <- c(-1 / Bbar, Abar / Bbar^2)
  S2 <- stats::cov(cbind(A, B))
  micro_se <- sqrt(drop(t(g2) %*% S2 %*% g2) / n)

  # gap = micro - macro = 1 - Abar/Bbar - rbar; gradient adds -1 for rbar
  g3 <- c(g2, -1)
  S3 <- stats::cov(cbind(A, B, rv))
  gap_se <- sqrt(drop(t(g3) %*% S3 %*% g3) / n)

  list(r = r, macro = macro, macro_se = macro_se,
       micro = micro, micro_se = micro_se, gap_se = gap_se, n = n)
}

#' Adjacent held-out gains and the epsilon-adequacy selection
#'
#' Paired comparison of held-out fit across adjacent topic counts, following
#' Proposition 2 (iv) and Definition 1 of Lewis and Grossetti (2026). For
#' each grid point \eqn{K} with successor \eqn{succ(K)} (the next point of
#' the actual grid), the per-document gain is
#' \eqn{\Delta_j(K) = R^{2,ho}_j(succ(K)) - R^{2,ho}_j(K)}{Delta_j(K) = R2_j(succ(K)) - R2_j(K)},
#' computed on the same evaluation documents. The table reports the mean
#' gain, its standard deviation, the one-sided upper confidence bound
#' \eqn{\bar\Delta + z_{1-\alpha}\,\hat\sigma_\Delta/\sqrt{J_+}}{mean + z * sd/sqrt(J+)},
#' and the paired test of no improvement. The selection
#' \eqn{\hat K_{D,\varepsilon,\alpha}}{K_hat} is the smallest \eqn{K} whose
#' upper bound is at most \eqn{\varepsilon}{epsilon}.
#'
#' @param ho An [optop_index_holdout()] result.
#' @param metric One of the metrics evaluated in `ho` (default
#'   `"deviance"`).
#' @param epsilon Positive tolerance of the adequacy rule (default `0.01`).
#' @param alpha One-sided level of the upper bound and the improvement test
#'   (default `0.05`).
#'
#' @details
#' The rule is a stopping rule over dependent one-sided tests: adjacent
#' gains share evaluation documents and overlapping models, so `alpha`
#' applies to each comparison separately, not to the selection event.
#' Inference reported after selection is subject to the usual post-selection
#' caveats; a conservative practice is to select \eqn{K} on one held-out
#' sample and report fit indices and moment diagnostics on a second,
#' disjoint evaluation sample. The rule is a practical sufficiency
#' criterion, not an estimator of a latent true number of topics.
#'
#' @return A list with
#' - `gains`: `data.table` with one row per adjacent pair: `K`, `succ_K`,
#'   `n` (paired documents), `gain` (mean), `sd`, `upper` (one-sided bound),
#'   `z`, and `pval` (one-sided, improvement).
#' - `k_hat`: the selected topic count (`NA` when no grid point qualifies).
#' - `metric`, `epsilon`, `alpha`.
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2026). Goodness-of-fit indices and
#' diagnostics for topic models. Working paper.
#'
#' @seealso [optop_index_holdout()]
#'
#' @export
optop_gain_table <- function(ho, metric = "deviance", epsilon = 0.01,
                             alpha = 0.05) {
  if (!inherits(ho, "optop_holdout")) {
    stop("ho must be an optop_index_holdout() result")
  }
  if (!metric %in% ho$metrics) {
    stop("metric was not evaluated in this holdout run")
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("epsilon must be a single positive number")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0, 1)")
  }

  S <- ho$scores[[metric]]
  K <- ho$K
  if (length(K) < 2L) stop("the gain table needs at least two topic counts")
  z1 <- stats::qnorm(1 - alpha)

  rows <- vector("list", length(K) - 1L)
  for (i in seq_len(length(K) - 1L)) {
    # Delta_j(K) = R2_j(succ(K)) - R2_j(K) on the same documents
    delta <- S[, i + 1L] - S[, i]
    delta <- delta[!is.na(delta)]
    n <- length(delta)
    m <- mean(delta)
    s <- stats::sd(delta)
    zstat <- sqrt(n) * m / s
    rows[[i]] <- data.table::data.table(
      K = K[i], succ_K = K[i + 1L], n = n,
      gain = m, sd = s,
      upper = m + z1 * s / sqrt(n),
      z = zstat,
      pval = stats::pnorm(zstat, lower.tail = FALSE)
    )
  }
  gains <- data.table::rbindlist(rows)

  ok <- gains$upper <= epsilon
  k_hat <- if (any(ok)) gains$K[which(ok)[1L]] else NA_real_

  list(gains = gains, k_hat = k_hat,
       metric = metric, epsilon = epsilon, alpha = alpha)
}
