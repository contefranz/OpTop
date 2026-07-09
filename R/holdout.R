#' Held-out goodness-of-fit indices with conditional inference
#'
#' Answer the question: does the chosen topic model also describe documents
#' it has never seen? Indices computed on the estimation corpus are
#' descriptive and optimistic, because the same documents shaped the topics
#' being evaluated. This function evaluates the same goodness-of-fit
#' indices on an independent *evaluation* sample instead, following
#' Section 3.7 of Lewis and Grossetti (2026), and attaches the uncertainty
#' measures that make held-out fit values comparable across topic counts:
#' a confidence interval for the average document's fit and standard errors
#' for the pooled index and the Micro-Macro gap.
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
#' **How it works.** The workflow has five steps, all handled internally:
#' 1. *Split upstream.* You fit the model grid on a training corpus and keep
#'    a separate set of evaluation documents that played no role in the
#'    fitting; the training baseline comes from
#'    `optop_make_baseline(dtm_train)`.
#' 2. *Fold in.* For every model and every evaluation document, the
#'    document's topic weights are estimated while the trained topics stay
#'    fixed (the fold-in step; see the engine notes below). Combined with
#'    the trained topics, the weights give the fitted word probabilities of
#'    each evaluation document.
#' 3. *Rebuild the common support.* The harmonized rare-word partition is
#'    reconstructed on the evaluation corpus, with the union running over
#'    the whole model grid plus the training baseline, exactly as in the
#'    in-sample analysis, so every topic count is scored on the same bins.
#' 4. *Score.* Each evaluation document contributes one fitted and one
#'    baseline discrepancy per metric, and its document-level index
#'    \eqn{1 - D_j(K)/D_j(\mathrm{null})}{1 - D_j(K)/D_j(null)}: the share
#'    of the no-topics discrepancy that the model removes on that document.
#' 5. *Aggregate.* The Macro index averages the document scores; the Micro
#'    index pools the discrepancies first. Because the evaluation documents
#'    are independent of one another once the training fits are fixed,
#'    both come with valid standard errors.
#'
#' **Held-out target.** The scores implement held-out-document
#' reconstruction: the same evaluation document is used to infer its topic
#' mixture and to score its fit, so the model is judged on how well it can
#' rebuild a document after reading it once. For the stricter
#' held-out-token completion target, split each evaluation document upstream
#' and pass the fold-in tokens through the engine directly.
#'
#' **Out-of-support convention.** Evaluation documents may contain words the
#' models were never trained on. Those tokens cannot be scored word by word:
#' they are assigned to the observed min-bin of their document (through a
#' synthetic always-rare column with zero fitted and baseline probability),
#' and document lengths \eqn{L_j} keep every token. Alternatively, smooth
#' the training baseline (`optop_make_baseline(dtm_train, smooth_lambda)`)
#' and report sensitivity.
#'
#' **Inference.** For each \eqn{K}, the Macro index is the mean of the
#' document scores over the non-degenerate evaluation documents
#' \eqn{J_+}{J+}; its standard error uses the \eqn{1/(J_+ - 1)}{1/(J+ - 1)}
#' sample variance and the interval is
#' \eqn{\mathrm{Macro} \pm z_{1-\alpha/2}\, \hat\sigma / \sqrt{J_+}}{Macro +/- z * sigma_hat / sqrt(J+)}
#' (Proposition 2 of the paper). Read the interval as the range of average
#' held-out fits compatible with the evaluation sample, for this training
#' fit: with a different evaluation sample from the same population, the
#' interval would cover the true average about `conf` of the time. The
#' Micro index and the Micro-Macro gap are smooth functions of conditional
#' means, and their standard errors follow by the delta method. All
#' inference is conditional on the training sample: the target is the
#' expected per-document fit under the training-fitted global objects, not
#' an unconditional population quantity.
#'
#' **Documents dropped.** Evaluation documents with zero aligned tokens
#' cannot be folded in and are dropped with a warning.
#'
#' @return An object of class `optop_holdout`: a list with
#' - `summary`: `data.table` with one row per (metric, K). Reading the
#'   columns: `macro` is the average document's held-out fit; `macro_se`
#'   its standard error; `ci_lo`, `ci_hi` the bounds of the `conf`-level
#'   interval for the average fit; `micro` the pooled index, in which
#'   documents with a large baseline discrepancy weigh more; `micro_se` its
#'   delta-method standard error; `gap` equals `micro - macro`, positive
#'   when fit concentrates in high-discrepancy (often long or atypical)
#'   documents; `gap_se` its standard error; `n_docs` the number of
#'   non-degenerate evaluation documents behind every estimate.
#' - `scores`: named list (one element per metric) of numeric matrices,
#'   evaluation documents by \eqn{K}, holding the per-document held-out
#'   indices (`NA` on degenerate documents when `stabilize = 0`). These are
#'   the scores [optop_gain_table()] pairs across topic counts.
#' - `d_model`, `d_null`: same shape, the per-document fitted and baseline
#'   discrepancies, for fit-versus-baseline diagnostics.
#' - `K`, `metrics`, `conf`, `stabilize`, `c`: the evaluation design.
#'
#' @examples
#' \donttest{
#' # a small synthetic corpus with 6 known topics
#' rdir <- function(n, k, a) {
#'   g <- matrix(stats::rgamma(n * k, shape = a), nrow = n)
#'   g / rowSums(g)
#' }
#' set.seed(42)
#' corpus <- sim_dfm(DTW = rdir(120, 6, 0.4), TWW = rdir(6, 300, 0.1),
#'                   doc_length = rep(400, 120), seed = 1)
#'
#' # split: 90 training documents, 30 evaluation documents
#' dtm_train <- corpus[1:90, ]
#' dtm_eval <- corpus[91:120, ]
#'
#' # fit a small grid on the TRAINING split only
#' models <- lapply(c(4, 6, 8), function(k) {
#'   topicmodels::LDA(dtm_train, k = k, control = list(seed = 100 + k))
#' })
#'
#' # held-out evaluation with the TRAINING baseline
#' base_tr <- optop_make_baseline(dtm_train)
#' ho <- optop_index_holdout(models, dtm_eval, base_tr)
#' ho$summary
#'
#' # when do extra topics stop paying for themselves?
#' gt <- optop_gain_table(ho, epsilon = 0.01)
#' gt$gains
#' gt$k_hat
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
#' Answer the question: when do additional topics stop paying for
#' themselves? Comparing two topic counts on held-out fit is a paired
#' problem: the same evaluation documents are scored under both models, so
#' the informative quantity is the per-document *gain*, the difference of
#' the two held-out scores on the same document. This function computes the
#' paired gains between adjacent points of the topic grid, attaches their
#' uncertainty, and applies the paper's selection rule: stop at the
#' smallest \eqn{K} beyond which the certified improvement falls below a
#' tolerance you chose in advance (Proposition 2 (iv) and Definition 1 of
#' Lewis and Grossetti 2026).
#'
#' @param ho An [optop_index_holdout()] result.
#' @param metric One of the metrics evaluated in `ho` (default
#'   `"deviance"`, the paper's primary index).
#' @param epsilon Positive tolerance of the adequacy rule: the held-out
#'   improvement you consider too small to justify an extra topic, on the
#'   proportional scale of the index (default `0.01`, one percentage point
#'   of discrepancy reduction). Fix it before looking at the results.
#' @param alpha One-sided level of the upper bound and the improvement test
#'   (default `0.05`).
#'
#' @details
#' **How it works.**
#' 1. For each grid point \eqn{K} with successor \eqn{succ(K)} (the next
#'    point of the actual grid, not mechanically \eqn{K + 1}), the
#'    per-document gain is
#'    \eqn{\Delta_j(K) = R^{2,ho}_j(succ(K)) - R^{2,ho}_j(K)}{Delta_j(K) = R2_j(succ(K)) - R2_j(K)},
#'    computed on the same evaluation documents.
#' 2. The mean gain estimates how much held-out fit the extra topics buy;
#'    its one-sided upper confidence bound
#'    \eqn{\bar\Delta + z_{1-\alpha}\,\hat\sigma_\Delta/\sqrt{J_+}}{mean + z * sd/sqrt(J+)}
#'    says: with confidence \eqn{1-\alpha}{1 - alpha}, the true gain is no
#'    larger than this.
#' 3. The selection \eqn{\hat K_{D,\varepsilon,\alpha}}{k_hat} is the
#'    smallest \eqn{K} whose upper bound is at most
#'    \eqn{\varepsilon}{epsilon}: even in the optimistic reading of the
#'    data, moving past that \eqn{K} buys less than the tolerance.
#'
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
#' - `gains`: `data.table` with one row per adjacent pair. Reading the
#'   columns: `K` and `succ_K` are the pair compared; `n` the number of
#'   paired documents; `gain` the mean per-document improvement from `K` to
#'   `succ_K`; `sd` its standard deviation across documents; `upper` the
#'   one-sided upper confidence bound of the mean gain (the quantity the
#'   selection rule thresholds); `z` and `pval` the paired test of no
#'   improvement (small `pval` = the extra topics helped on average).
#' - `k_hat`: the selected topic count, the smallest `K` with
#'   `upper <= epsilon` (`NA` when no grid point qualifies, meaning the
#'   grid ended before the gains died out).
#' - `metric`, `epsilon`, `alpha`: the design.
#'
#' @examples
#' \donttest{
#' # continue from the optop_index_holdout() example
#' rdir <- function(n, k, a) {
#'   g <- matrix(stats::rgamma(n * k, shape = a), nrow = n)
#'   g / rowSums(g)
#' }
#' set.seed(42)
#' corpus <- sim_dfm(DTW = rdir(120, 6, 0.4), TWW = rdir(6, 300, 0.1),
#'                   doc_length = rep(400, 120), seed = 1)
#' dtm_train <- corpus[1:90, ]
#' dtm_eval <- corpus[91:120, ]
#' models <- lapply(c(4, 6, 8), function(k) {
#'   topicmodels::LDA(dtm_train, k = k, control = list(seed = 100 + k))
#' })
#' ho <- optop_index_holdout(models, dtm_eval,
#'                           optop_make_baseline(dtm_train))
#'
#' gt <- optop_gain_table(ho, epsilon = 0.01, alpha = 0.05)
#' gt$gains   # certified improvement per grid step
#' gt$k_hat   # smallest K after which gains fall below epsilon
#' }
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2026). Goodness-of-fit indices and
#' diagnostics for topic models. Working paper.
#'
#' @seealso [optop_index_holdout()], [optop_moment_test()]
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
