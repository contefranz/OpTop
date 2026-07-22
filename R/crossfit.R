# V-fold cross-fitted discrepancy indices (paper Appendix B): partition the
# documents into V folds, refit the model grid on each fold's complement,
# and score every fold through the held-out machinery of R/holdout.R, so
# each document is evaluated exactly once, by models that never saw it.
# Point estimates pool the per-document out-of-fold scores; standard errors
# come from the V fold replicates, which carry the fold-to-fold training
# variation a single split cannot see. The per-fold work reuses
# optop_index_holdout() wholesale: the training baseline, the out-of-support
# convention, fold-in, the held-out partition, and the floor all behave
# exactly as documented there.

#' V-fold cross-fitted discrepancy indices
#'
#' Cross-fit the held-out discrepancy indices of [optop_index_holdout()]:
#' split the corpus into `V` folds, refit the model grid on each fold's
#' complement with a user-supplied fitting function, score the held-out
#' fold, and pool. Every document is scored exactly once, by models
#' trained without it, so the pooled indices are out-of-sample for the
#' whole corpus rather than for one evaluation split; this is the
#' implementation the paper recommends in practice (Appendix B).
#'
#' @section The fitting function:
#' `fit_fun(dtm_train, k)` receives the training complement of the fold
#' (rows of `dtm` in the class you passed) and one topic count, and must
#' return a fitted model OpTop supports: anything
#' [optop_as_theta_phi()] adapts *and* [optop_fold_in()] can project, so
#' fits from topicmodels, seededlda, NLPstudio, and text2vec via
#' [optop_warplda()] all qualify. A bare [optop_model()] does not: it has
#' no engine to fold unseen documents in. Seeding the engine is
#' `fit_fun`'s own business; make it depend on `k` (and, if you want
#' fold-level reproducibility, on the training row names) for
#' reproducible fits.
#'
#' @section Pooling and standard errors:
#' Point estimates are computed at the document level on the pooled
#' out-of-fold discrepancies, under the same null-discrepancy floor and
#' stabilization conventions as [optop_index_holdout()]. Standard errors
#' are **fold-replicate (cluster) estimates**: the V per-fold values of
#' each quantity (Macro, Micro, the Micro-Macro gap, and the adjacent
#' gains) are treated as exchangeable replicates and the standard error
#' of their mean is `sd / sqrt(V)`. This estimator sees the fold-to-fold
#' variation induced by refitting, which the per-document standard error
#' of a single split ignores; it is the conservative choice the paper's
#' appendix recommends, at the price of only V replicates (a warning is
#' issued when folds are small). The epsilon-adequacy selection of
#' [optop_gain_table()] is applied to the cross-fitted gains with these
#' fold-replicate bounds; the result is reported per metric and
#' recomputable at other tolerances through [plot.optop_crossfit()].
#'
#' @param dtm A document-term matrix of raw counts (a `quanteda::dfm`, a
#'   `dgCMatrix`, or a base matrix with dimnames). Sharded
#'   [optop_corpus()] objects are not supported: cross-fitting refits the
#'   grid `V` times, which is engine-scale work outside the sharded
#'   evaluation contract.
#' @param K Integer vector of topic counts, strictly increasing, all at
#'   least 2. The grid is refit in every fold.
#' @param fit_fun A function `(dtm_train, k)` returning a fitted model;
#'   see the dedicated section.
#' @param V Number of folds, between 2 and `floor(J / 2)`; the default 5
#'   follows the paper's appendix.
#' @param metrics Character subset of `c("se", "chisq", "deviance")`.
#' @param c Positive rare-word threshold constant of the harmonized
#'   partition, as in [optop_make_partition()].
#' @param conf Confidence level of the reported intervals.
#' @param stabilize Nonnegative stabilization constant kappa, as in
#'   [optop_index_holdout()].
#' @param min_null Null-discrepancy floor; `NULL` (the default) resolves
#'   to `c`, and `0` restores the strict-positivity rule.
#' @param stratify Logical; when `TRUE` (the default) folds are balanced
#'   on document length (documents are ordered by length and each
#'   consecutive block of `V` is dealt across folds at random), so no
#'   fold concentrates the short documents. When `FALSE` the assignment
#'   is simple random.
#' @param seed Optional integer seeding the fold assignment (the only
#'   randomness this function itself draws).
#' @param n_threads OpenMP threads for the compiled scoring kernels;
#'   never changes results.
#' @param verbose Logical; when `TRUE`, progress alerts per fold.
#' @param ... Passed to [optop_fold_in()] through the per-fold scoring
#'   (engine-specific fold-in controls).
#'
#' @return An object of class `optop_crossfit`: a list with
#'   \describe{
#'     \item{`summary`}{One row per metric and K: pooled `macro` with
#'       fold-replicate `macro_se` and its interval, pooled `micro` and
#'       `micro_se`, the Micro-Macro `gap` and `gap_se`, `n_docs` in the
#'       pooled J+, and `n_null_excluded`.}
#'     \item{`gains`}{Per metric, the cross-fitted adjacent gains with
#'       fold-replicate upper bounds and the selected `k_hat` at the
#'       default tolerance (`epsilon = 0.01`, `alpha = 0.05`).}
#'     \item{`scores`}{Per metric, the J x length(K) matrix of
#'       out-of-fold document indices (NA for floor-excluded documents).}
#'     \item{`d_model`, `d_null`}{The pooled out-of-fold discrepancies.}
#'     \item{`folds`}{The fold assignment, named by document.}
#'     \item{`fold_summary`}{The per-fold holdout summaries with a
#'       `fold` column, the replicates behind every standard error.}
#'     \item{`K`, `V`, `metrics`, `conf`, `stabilize`, `c`, `min_null`,
#'       `seed`}{The configuration.}
#'   }
#'
#' @examples
#' \donttest{
#' # simulate a corpus with known truth
#' rdirich <- function(n, k) {
#'   g <- matrix(stats::rgamma(n * k, shape = 1), n, k)
#'   g / rowSums(g)
#' }
#' theta <- rdirich(60, 3)
#' phi <- rdirich(3, 120)
#' colnames(phi) <- sprintf("w%03d", 1:120)
#' dtm <- sim_dfm(theta, phi, doc_length = 150, seed = 7)
#'
#' fit_vem <- function(dtm_train, k) {
#'   topicmodels::LDA(quanteda::convert(dtm_train, to = "topicmodels"),
#'                    k = k, method = "VEM",
#'                    control = list(seed = 1000 + k))
#' }
#' cf <- optop_crossfit(dtm, K = 2:4, fit_fun = fit_vem, V = 3,
#'                      metrics = "deviance", seed = 42)
#' cf
#' plot(cf)
#' }
#'
#' @seealso [optop_index_holdout()] for the single-split machinery this
#'   builds on, [optop_gain_table()] for the selection rule applied to
#'   the cross-fitted gains, [plot.optop_crossfit()].
#'
#' @export
optop_crossfit <- function(dtm, K, fit_fun, V = 5,
                           metrics = c("se", "chisq", "deviance"),
                           c = 1, conf = 0.95, stabilize = 0,
                           min_null = NULL, stratify = TRUE, seed = NULL,
                           n_threads = 1L, verbose = FALSE, ...) {
  if (.optop_is_corpus(dtm)) {
    stop(paste("optop_crossfit() does not accept a sharded optop_corpus():",
               "cross-fitting refits the model grid V times, which is",
               "engine-scale work; run it on a corpus that fits one matrix,",
               "or use optop_index_holdout() with a fixed training split."))
  }
  if (is.null(rownames(dtm))) {
    stop("dtm needs document names; the pooled scores are aligned by name")
  }
  J <- nrow(dtm)
  K <- as.integer(K)
  if (anyNA(K) || any(K < 2) || is.unsorted(K, strictly = TRUE)) {
    stop("K must be a strictly increasing integer vector with all K >= 2")
  }
  if (!is.function(fit_fun) || length(formals(fit_fun)) < 2L) {
    stop("fit_fun must be a function of (dtm_train, k)")
  }
  V <- as.integer(V)
  if (length(V) != 1L || is.na(V) || V < 2L || V > J %/% 2L) {
    stop("V must be a single integer between 2 and floor(J / 2)")
  }
  metrics <- intersect(c("se", "chisq", "deviance"), metrics)
  if (!length(metrics)) stop("no valid metric requested")
  min_null <- .optop_resolve_min_null(min_null, c)
  if (J %/% V < 30L) {
    warning("fewer than 30 documents per fold: the fold Macro means are ",
            "noisy replicates and the fold-replicate standard errors ",
            "will be unstable", call. = FALSE)
  }

  L <- as.numeric(Matrix::rowSums(dtm))
  if (!is.null(seed)) set.seed(seed)
  folds <- .optop_crossfit_folds(L, V, stratify)
  names(folds) <- rownames(dtm)

  if (verbose) {
    cli::cli_h2("V-fold cross-fitting")
    cli::cli_alert_info(paste0(
      "{J} documents, V = {V} folds, K grid {{", paste(K, collapse = ", "),
      "}}, metrics: ", paste(metrics, collapse = ", ")))
  }

  n_k <- length(K)
  scores <- d_model <- d_null <- stats::setNames(lapply(metrics, function(m) {
    matrix(NA_real_, J, n_k, dimnames = list(rownames(dtm), K))
  }), metrics)
  fold_rows <- vector("list", V)

  for (v in seq_len(V)) {
    idx_ev <- which(folds == v)
    dtm_tr <- dtm[-idx_ev, , drop = FALSE]
    dtm_ev <- dtm[idx_ev, , drop = FALSE]
    if (verbose) {
      cli::cli_alert_info(paste0(
        "fold {v}/{V}: fitting {n_k} model{?s} on {nrow(dtm_tr)} documents,",
        " scoring {nrow(dtm_ev)}"))
    }

    models_v <- lapply(K, function(k) {
      tryCatch(fit_fun(dtm_tr, k), error = function(e) {
        stop("fit_fun failed for K = ", k, " in fold ", v, ": ",
             conditionMessage(e), call. = FALSE)
      })
    })
    baseline_v <- optop_make_baseline(dtm_tr)
    ho_v <- suppressMessages(
      optop_index_holdout(models_v, dtm_ev, baseline_v, c = c,
                          metrics = metrics, conf = conf,
                          stabilize = stabilize, n_threads = n_threads,
                          min_null = min_null, ...)
    )
    if (!identical(as.integer(ho_v$K), K)) {
      stop("fit_fun returned models with topic counts {",
           paste(ho_v$K, collapse = ", "), "} where K = {",
           paste(K, collapse = ", "), "} was requested (fold ", v, ")")
    }

    for (m in metrics) {
      rows_ev <- rownames(ho_v$scores[[m]])
      scores[[m]][rows_ev, ] <- ho_v$scores[[m]]
      d_model[[m]][rows_ev, ] <- ho_v$d_model[[m]]
      d_null[[m]][rows_ev, ] <- ho_v$d_null[[m]]
    }
    fold_rows[[v]] <- data.table::data.table(fold = v, ho_v$summary)
  }
  fold_summary <- data.table::rbindlist(fold_rows)

  # pooled point estimates at the document level, fold-replicate variances
  z <- stats::qnorm(1 - (1 - conf) / 2)
  rows <- list()
  for (m in metrics) {
    for (i in seq_len(n_k)) {
      st <- .optop_holdout_stats(d_model[[m]][, i], d_null[[m]][, i],
                                 stabilize, z, min_null)
      reps <- fold_summary[fold_summary$metric == m &
                             fold_summary$K == K[i], ]
      se_of <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2L) return(NA_real_)
        stats::sd(x) / sqrt(length(x))
      }
      macro_se <- se_of(reps$macro)
      rows[[length(rows) + 1L]] <- data.table::data.table(
        metric = m, K = K[i],
        macro = st$macro, macro_se = macro_se,
        ci_lo = st$macro - z * macro_se,
        ci_hi = st$macro + z * macro_se,
        micro = st$micro, micro_se = se_of(reps$micro),
        gap = st$micro - st$macro, gap_se = se_of(reps$gap),
        n_docs = st$n, n_null_excluded = st$n_excluded
      )
    }
  }
  summary <- data.table::rbindlist(rows)

  # the pooled floor report, once per metric
  for (m in metrics) {
    n_excl <- max(summary$n_null_excluded[summary$metric == m])
    .optop_report_null_floor(n_excl, n_excl / J, min_null, m)
  }

  out <- list(summary = summary, scores = scores,
              d_model = d_model, d_null = d_null,
              folds = folds, fold_summary = fold_summary,
              K = K, V = V, metrics = metrics, conf = conf,
              stabilize = stabilize, c = c, min_null = min_null,
              seed = seed)
  out$gains <- stats::setNames(lapply(metrics, function(m) {
    .optop_crossfit_gains(out, m, epsilon = 0.01, alpha = 0.05)
  }), metrics)
  class(out) <- c("optop_crossfit", "list")
  out
}

#' Length-balanced fold assignment
#'
#' Orders the documents by length and deals each consecutive block of V
#' documents across the V folds at random, so fold sizes differ by at
#' most one and every fold spans the length distribution. With
#' `stratify = FALSE` the assignment is a simple random permutation of a
#' balanced label vector.
#'
#' @param L Numeric vector of document lengths.
#' @param V Number of folds.
#' @param stratify Logical.
#'
#' @return Integer fold labels in `1:V`, one per document.
#'
#' @keywords internal
.optop_crossfit_folds <- function(L, V, stratify) {
  J <- length(L)
  folds <- integer(J)
  if (stratify) {
    ord <- order(L, seq_len(J))
    pos <- 1L
    while (pos <= J) {
      block <- ord[pos:min(pos + V - 1L, J)]
      folds[block] <- sample(V, length(block))
      pos <- pos + V
    }
  } else {
    folds <- sample(rep_len(seq_len(V), J))
  }
  folds
}

#' Cross-fitted adjacent gains with fold-replicate bounds
#'
#' The cross-fitted analogue of [optop_gain_table()]: per-document
#' adjacent gains are pooled over all out-of-fold scores, and the
#' one-sided upper bound uses the fold-replicate standard error (the V
#' per-fold mean gains as replicates) instead of the per-document one.
#' The epsilon-adequacy rule then selects the first K whose upper bound
#' falls below `epsilon`.
#'
#' @param x An `optop_crossfit` result.
#' @param metric One evaluated metric.
#' @param epsilon,alpha The adequacy tolerance and one-sided level.
#'
#' @return A list with the `gains` table (`K`, `succ_K`, `n`, `gain`,
#'   `se_fold`, `upper`), the selected `k_hat` (NA when no step
#'   qualifies), and the configuration.
#'
#' @keywords internal
.optop_crossfit_gains <- function(x, metric, epsilon, alpha) {
  S <- x$scores[[metric]]
  K <- x$K
  if (length(K) < 2L) {
    return(list(gains = data.table::data.table(), k_hat = NA_real_,
                metric = metric, epsilon = epsilon, alpha = alpha))
  }
  z1 <- stats::qnorm(1 - alpha)
  rows <- vector("list", length(K) - 1L)
  for (i in seq_len(length(K) - 1L)) {
    delta <- S[, i + 1L] - S[, i]
    ok <- !is.na(delta)
    m <- mean(delta[ok])
    gv <- tapply(delta[ok], x$folds[ok], mean)
    se_fold <- if (length(gv) < 2L) NA_real_ else {
      stats::sd(gv) / sqrt(length(gv))
    }
    rows[[i]] <- data.table::data.table(
      K = K[i], succ_K = K[i + 1L], n = sum(ok),
      gain = m, se_fold = se_fold,
      upper = m + z1 * se_fold
    )
  }
  gains <- data.table::rbindlist(rows)
  ok_step <- !is.na(gains$upper) & gains$upper <= epsilon
  k_hat <- if (any(ok_step)) gains$K[which(ok_step)[1L]] else NA_real_
  list(gains = gains, k_hat = k_hat,
       metric = metric, epsilon = epsilon, alpha = alpha)
}