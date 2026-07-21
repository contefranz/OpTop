#' Held-out moment-based specification tests
#'
#' Answer the question: which parts of the vocabulary does the model still
#' get systematically wrong? A fit index says how much error remains, not
#' where. This function groups the vocabulary into researcher-specified,
#' labeled word sets (the *instruments*: for example frequency strata, or
#' strata of a training fit score) and tests whether the model consistently
#' over- or under-predicts any group on documents it has never seen,
#' following Section 4 of Lewis and Grossetti (2026).
#'
#' Formally, for each evaluation document the residual probability vector is
#' \eqn{\varepsilon_j = \mathbf d_j - \hat{\mathbf p}^{K,tr}_j}{e_j = d_j - p_j},
#' observed minus fitted word probability under the training-fitted model;
#' the moment vector projects it onto a training-built instrument matrix,
#' \eqn{\mathbf g_j = \mathbf Z\,\varepsilon_j}{g_j = Z e_j}, and the test
#' asks whether the moments average to zero across evaluation documents.
#' Positive components mean the model under-predicts that direction of the
#' vocabulary.
#'
#' @param models A list of topic models fitted on the \emph{training}
#'   corpus (any fold-in-capable engine; see [optop_index_holdout()]).
#' @param dtm_eval A counts document-term matrix of the evaluation
#'   documents; columns are matched to the training vocabulary by name.
#' @param dtm_train The counts document-term matrix of the \emph{training}
#'   corpus, used only to build the instruments (word frequencies and, for
#'   `type = "fit"`, the training word-level fit score); never touched by
#'   the evaluation residuals.
#' @param type Instrument family:
#'   * `"contrast"` (Test 1): one row contrasting the top against the bottom
#'     training-frequency quintile; the scalar `t` statistic applies.
#'   * `"strata"` (Test 2): `bins - 1` rows comparing each training-frequency
#'     stratum against the highest-frequency reference stratum.
#'   * `"fit"` (Test 3): `strata - 1` rows comparing strata of the training
#'     word-level deviance index against the best-fit reference stratum; the
#'     instrument matrix depends on \eqn{K}.
#' @param bins Number of frequency strata for `type = "strata"` (default 5).
#' @param strata Number of fit strata for `type = "fit"` (default 5).
#' @param min_doc_freq Minimum training document frequency for a word to
#'   enter the fit strata (default 5); excluded words carry zero instrument
#'   entries, which leaves every row sum at zero.
#' @param adjust Multiple-testing adjustment of the marginal per-stratum
#'   t tests: `"none"` (default), `"bonferroni"`, or `"BH"`.
#' @param baseline Optional \emph{training} baseline for the fit score
#'   (default `NULL` computes it from `dtm_train`).
#' @param n_threads Integer; OpenMP threads of the compiled kernels used for
#'   the training fit score.
#' @param ... Passed to the engine's fold-in routine.
#'
#' @details
#' **How it works.**
#' 1. *Split upstream.* Models are fitted on the training corpus;
#'    the evaluation documents stay out of the fitting entirely.
#' 2. *Build the instruments from training data only.* Each instrument row
#'    contrasts one word group against a reference group, so its moment
#'    reads as "average residual mass in this group minus the reference".
#'    Because the groups are fixed before any evaluation residual is seen,
#'    the test is not fishing in its own evaluation data.
#' 3. *Score every evaluation document.* Topic weights are folded in
#'    (holding the trained topics fixed), the residual vector is formed,
#'    and its projection on the instruments gives one small moment vector
#'    per document.
#' 4. *Test.* If the model has no systematic bias along the chosen
#'    directions, the moments should average to zero; the Wald statistic
#'    measures how far from zero the averages are relative to their
#'    sampling noise.
#'
#' **Construction.** All instruments are built from the training sample
#' only, before any evaluation residual is examined, and every row is
#' exactly mean-zero with \eqn{\|z\|_1 = 2}{||z||_1 = 2}. The residuals use
#' the evaluation counts on the training vocabulary with document lengths
#' restricted to the aligned tokens, so that observed and fitted vectors are
#' both probability vectors on the training support and
#' \eqn{\mathbf 1^\top \varepsilon_j = 0}{sum(e_j) = 0} holds identically;
#' out-of-support tokens play no role in the residual moments.
#'
#' **Statistics.** With `n` evaluation documents, the Wald statistic
#' \eqn{\mathcal W = n\,\bar{\mathbf g}^\top \hat\Sigma^{-1} \bar{\mathbf g}}{W = n * gbar' S^-1 gbar}
#' is asymptotically \eqn{\chi^2_q}{chi-square(q)} under the null of no
#' systematic residual bias (Proposition 3 of the paper), where
#' \eqn{\hat\Sigma}{S} is the sample covariance of the document moments. In
#' the scalar case the equivalent `t` statistic is standard normal and
#' \eqn{\mathcal W = t^2}{W = t^2}. Marginal per-stratum t tests localize a
#' rejection, with the requested adjustment.
#'
#' **Interpretation.** The null is a statement about the average held-out
#' residual balance of the fitted model, conditional on the training sample;
#' it is not equivalent to correct specification of the model class.
#' Estimation error and the reuse of the evaluation document in folding in
#' its topic weights generically make the moments small but nonzero, so with
#' many evaluation documents the test detects any residual imbalance,
#' whatever its source. Read rejections jointly with the magnitude of the
#' mean moments, which are effect sizes in probability-mass units. When many
#' instruments, topic counts, or families are examined, treat the p-values
#' as exploratory or adjust for multiple testing.
#'
#' @return An object of class `optop_moment_test`: a list with
#' - `summary`: `data.table` with one row per \eqn{K}. Reading the columns:
#'   `q` is the number of instrument rows (word-group contrasts); `wald`
#'   the joint test statistic with `df = q` and `pval` its p-value, where a
#'   small value indicates systematic residual bias along at least one
#'   contrast; `t` and `pval_t` are the equivalent scalar test when
#'   `q = 1` (`NA` otherwise).
#' - `moments`: per-\eqn{K} list with `g_bar` (the mean moments, effect
#'   sizes in probability-mass units: how much observed mass exceeds
#'   fitted mass along each contrast), `sigma` (their covariance), `n`, and
#'   `marginal` (`data.table` localizing a rejection: one row per stratum
#'   with `estimate`, `se`, `t`, `pval`, and the adjusted `pval_adj`).
#' - `instruments`: the instrument matrix and the stratum assignment of
#'   every vocabulary word (per \eqn{K} for `type = "fit"`).
#' - `type`, `adjust`, `K`: the design.
#'
#' @examples
#' \donttest{
#' # a small synthetic corpus, split into training and evaluation
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
#'
#' # does residual bias vary across frequency strata?
#' mt <- optop_moment_test(models, dtm_eval, dtm_train, type = "strata",
#'                         bins = 4)
#' mt$summary                  # one Wald test per K
#' mt$moments[["6"]]$marginal  # which stratum drives it, at K = 6
#' }
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2026). Goodness-of-fit indices and
#' diagnostics for topic models. Working paper.
#'
#' Hansen, L. P. (1982). Large sample properties of generalized method of
#' moments estimators. *Econometrica*, 50(4), 1029--1054.
#'
#' @seealso [optop_index_holdout()], [optop_gain_table()]
#'
#' @export
optop_moment_test <- function(models, dtm_eval, dtm_train,
                              type = c("contrast", "strata", "fit"),
                              bins = 5, strata = 5, min_doc_freq = 5,
                              adjust = c("none", "bonferroni", "BH"),
                              baseline = NULL, n_threads = 1L, ...) {
  type <- match.arg(type)
  adjust <- match.arg(adjust)
  stopifnot(length(models) >= 1)

  if (is.null(baseline)) baseline <- optop_make_baseline(dtm_train)

  # sharded evaluation corpora stream one shard at a time; the per-document
  # moment vectors row-bind exactly because every g_j depends on document j
  # alone once the training objects are fixed
  eval_shards <- if (.optop_is_corpus(dtm_eval)) {
    lapply(seq_len(dtm_eval$n_shards),
           function(s) .optop_corpus_shard(dtm_eval, s))
  } else {
    list(dtm_eval)
  }
  preps <- lapply(eval_shards, function(sh) {
    .optop_holdout_prepare(models, sh, baseline, ...)
  })
  prep <- preps[[1L]]
  vocab <- prep$vocab_tr
  W <- length(vocab)
  for (p in preps[-1L]) {
    if (!identical(p$K, prep$K)) {
      stop("the model grid changed across evaluation shards")
    }
  }

  if (.optop_is_corpus(dtm_train)) {
    if (!identical(.optop_corpus_vocab(dtm_train), vocab)) {
      stop(paste("dtm_train vocabulary/order differs from the models;",
                 "align the training dtm first"))
    }
    f_tr <- numeric(W)
    for (s in seq_len(dtm_train$n_shards)) {
      f_tr <- f_tr +
        as.numeric(Matrix::colSums(.optop_corpus_shard(dtm_train, s)))
    }
  } else {
    dtm_train <- methods::as(dtm_train, "CsparseMatrix")
    if (!identical(colnames(dtm_train), vocab)) {
      stop(paste("dtm_train vocabulary/order differs from the models;",
                 "align the training dtm first"))
    }
    f_tr <- as.numeric(Matrix::colSums(dtm_train))
  }

  # evaluation residual ingredients: d_j on the training support with
  # L_j restricted to aligned tokens, so sum_w d_jw = 1 and 1' e_j = 0
  n_ev <- sum(vapply(preps, function(p) nrow(p$dtm_aligned), numeric(1)))
  if (n_ev < 2L) stop("at least two evaluation documents are required")

  # training-frequency instruments are shared across the grid
  Z_fixed <- switch(type,
    contrast = .optop_z_contrast(f_tr, vocab),
    strata = .optop_z_strata(f_tr, vocab, bins),
    fit = NULL
  )

  K <- prep$K
  n_models <- length(prep$models)
  rows <- vector("list", n_models)
  moments <- stats::setNames(vector("list", n_models),
                             as.character(K))
  instruments <- if (type == "fit") {
    stats::setNames(vector("list", n_models), as.character(K))
  } else {
    Z_fixed
  }

  for (i_mod in seq_along(prep$models)) {
    Zi <- if (type == "fit") {
      .optop_z_fit(prep$models[[i_mod]], dtm_train, baseline, vocab,
                   strata, min_doc_freq, n_threads)
    } else {
      Z_fixed
    }
    if (type == "fit") instruments[[i_mod]] <- Zi

    # g_j = Z d_j - theta_j (phi Z^T): no J x W temporary is formed, and
    # shards contribute independent blocks of rows
    Zt <- Matrix::t(Zi$Z)
    G <- do.call(rbind, lapply(preps, function(p) {
      N_ev <- p$dtm_aligned
      L_ev <- Matrix::rowSums(N_ev)
      A <- as.matrix(N_ev %*% Zt) / L_ev
      tp <- p$tp_eval[[i_mod]]
      phi_tr <- tp$phi[, seq_len(W), drop = FALSE]
      B <- tp$theta %*% as.matrix(phi_tr %*% Zt)
      A - B
    }))

    q <- ncol(G)
    g_bar <- colMeans(G)
    sigma <- stats::cov(G)
    inv <- tryCatch(solve(sigma), error = function(e) {
      stop(paste("the moment covariance is numerically singular",
                 "(non-degeneracy assumption); reduce the number of",
                 "strata or enlarge the evaluation sample"))
    })
    wald <- n_ev * drop(t(g_bar) %*% inv %*% g_bar)
    pval <- stats::pchisq(wald, df = q, lower.tail = FALSE)

    se_b <- sqrt(diag(sigma) / n_ev)
    t_b <- g_bar / se_b
    p_b <- 2 * stats::pnorm(abs(t_b), lower.tail = FALSE)
    marginal <- data.table::data.table(
      stratum = Zi$labels,
      estimate = g_bar, se = se_b, t = t_b, pval = p_b,
      pval_adj = stats::p.adjust(p_b, method = if (adjust == "none") "none"
                                 else adjust)
    )

    t_stat <- if (q == 1L) sqrt(n_ev) * g_bar[1L] /
      stats::sd(G[, 1L]) else NA_real_
    rows[[i_mod]] <- data.table::data.table(
      K = K[i_mod], q = q, wald = wald, df = q, pval = pval,
      t = t_stat,
      pval_t = if (q == 1L) 2 * stats::pnorm(abs(t_stat),
                                             lower.tail = FALSE)
               else NA_real_
    )
    moments[[i_mod]] <- list(g_bar = g_bar, sigma = sigma, n = n_ev,
                             marginal = marginal)
  }

  out <- list(summary = data.table::rbindlist(rows),
              moments = moments, instruments = instruments,
              type = type, adjust = adjust, K = K)
  class(out) <- c("optop_moment_test", "list")
  out
}

#' Frequency-contrast instrument (Test 1)
#'
#' One mean-zero row over the vocabulary:
#' `z_w = 1(w in top quintile)/n_hi - 1(w in bottom quintile)/n_lo`, with the
#' quintiles taken over the ranked training frequencies (deterministic tie
#' handling by vocabulary position), so sum_w z_w = 0 and ||z||_1 = 2.
#'
#' @param f Training corpus frequencies, one per word.
#' @param vocab The vocabulary (column labels of the instrument).
#'
#' @return A list with the sparse `Z` (1 x W), stratum `labels`, and the
#'   per-word `assignment` vector.
#'
#' @keywords internal
.optop_z_contrast <- function(f, vocab) {
  W <- length(f)
  n_q <- max(1L, floor(W / 5))
  o <- order(f, seq_along(f))
  lo <- o[seq_len(n_q)]
  hi <- o[seq.int(W - n_q + 1L, W)]
  z <- numeric(W)
  z[hi] <- 1 / length(hi)
  z[lo] <- -1 / length(lo)
  assignment <- rep(NA_character_, W)
  assignment[lo] <- "low"
  assignment[hi] <- "high"
  Z <- Matrix::Matrix(matrix(z, nrow = 1, dimnames = list("hi_vs_lo", vocab)),
                      sparse = TRUE)
  list(Z = Z, labels = "hi_vs_lo", assignment = assignment)
}

#' Frequency-strata instruments (Test 2)
#'
#' bins - 1 mean-zero rows: row b compares the average residual mass of
#' frequency stratum b against the highest-frequency reference stratum B,
#' `Z_bw = 1(w in V_b)/|V_b| - 1(w in V_B)/|V_B|`. Strata are contiguous
#' chunks of the ranked training frequencies (equal sizes up to rounding,
#' deterministic tie handling).
#'
#' @param f Training corpus frequencies.
#' @param vocab The vocabulary.
#' @param bins Number of strata (B >= 2).
#'
#' @return As in `.optop_z_contrast()`, with `Z` of dimension (B-1) x W.
#'
#' @keywords internal
.optop_z_strata <- function(f, vocab, bins) {
  if (!is.numeric(bins) || length(bins) != 1L || bins < 2) {
    stop("bins must be a single integer >= 2")
  }
  bins <- as.integer(bins)
  W <- length(f)
  if (W < bins) stop("fewer words than frequency strata")
  o <- order(f, seq_along(f))
  cuts <- floor(seq(0L, W, length.out = bins + 1L))
  assignment <- integer(W)
  for (b in seq_len(bins)) {
    assignment[o[seq.int(cuts[b] + 1L, cuts[b + 1L])]] <- b
  }
  ref <- assignment == bins
  n_ref <- sum(ref)
  rows <- lapply(seq_len(bins - 1L), function(b) {
    z <- numeric(W)
    z[assignment == b] <- 1 / sum(assignment == b)
    z[ref] <- -1 / n_ref
    z
  })
  Z <- Matrix::Matrix(do.call(rbind, rows), sparse = TRUE)
  labels <- sprintf("freq_%d_vs_%d", seq_len(bins - 1L), bins)
  dimnames(Z) <- list(labels, vocab)
  list(Z = Z, labels = labels, assignment = as.character(assignment))
}

#' Fit-stratified instruments (Test 3)
#'
#' strata - 1 mean-zero rows comparing strata of the training word-level
#' deviance index against the best-fit reference stratum. Words with zero
#' training baseline deviance or training document frequency below
#' min_doc_freq are excluded from the strata and carry zero instrument
#' entries, which leaves every row sum at zero.
#'
#' @param model The training fit of the K under test (the score is
#'   K-specific).
#' @param dtm_train Training counts.
#' @param baseline Training baseline.
#' @param vocab The vocabulary.
#' @param strata Number of strata (S >= 2).
#' @param min_doc_freq Minimum training document frequency.
#' @param n_threads OpenMP threads for the engine sweep.
#'
#' @return As in `.optop_z_contrast()`, with `Z` of dimension (S-1) x W.
#'
#' @keywords internal
.optop_z_fit <- function(model, dtm_train, baseline, vocab, strata,
                         min_doc_freq, n_threads) {
  if (!is.numeric(strata) || length(strata) != 1L || strata < 2) {
    stop("strata must be a single integer >= 2")
  }
  strata <- as.integer(strata)
  W <- length(vocab)

  # training word-level deviance index s_K(w), unbinned support: the word
  # kernel needs only the document lengths from the partition object
  tp <- optop_as_theta_phi(.optop_materialize_model(model))
  pi_row <- baseline$pi_glob[vocab]
  if (.optop_is_corpus(dtm_train)) {
    L_tr <- unlist(lapply(seq_len(dtm_train$n_shards), function(s) {
      Matrix::rowSums(.optop_corpus_shard(dtm_train, s))
    }))
    docfreq <- numeric(W)
    for (s in seq_len(dtm_train$n_shards)) {
      docfreq <- docfreq +
        as.numeric(Matrix::colSums(.optop_corpus_shard(dtm_train, s) > 0))
    }
  } else {
    L_tr <- Matrix::rowSums(dtm_train)
    docfreq <- as.numeric(Matrix::colSums(dtm_train > 0))
  }
  part_stub <- list(L = L_tr)
  eng <- .optop_index_engine(tp$theta, tp$phi, dtm_train, part_stub, pi_row,
                             "deviance", "word", NULL, n_threads,
                             do_model = TRUE, do_null = TRUE)
  d_w <- eng$deviance$model
  d_null <- eng$deviance$null
  include <- d_null > 0 & docfreq >= min_doc_freq
  n_inc <- sum(include)
  if (n_inc < strata) {
    stop("fewer supported words than fit strata; lower min_doc_freq")
  }
  score <- 1 - d_w[include] / d_null[include]

  idx_inc <- which(include)
  o <- idx_inc[order(score, seq_along(score))]
  cuts <- floor(seq(0L, n_inc, length.out = strata + 1L))
  assignment <- rep(NA_integer_, W)
  for (s in seq_len(strata)) {
    assignment[o[seq.int(cuts[s] + 1L, cuts[s + 1L])]] <- s
  }
  ref <- !is.na(assignment) & assignment == strata
  n_ref <- sum(ref)
  rows <- lapply(seq_len(strata - 1L), function(s) {
    z <- numeric(W)
    sel <- !is.na(assignment) & assignment == s
    z[sel] <- 1 / sum(sel)
    z[ref] <- -1 / n_ref
    z
  })
  Z <- Matrix::Matrix(do.call(rbind, rows), sparse = TRUE)
  labels <- sprintf("fit_%d_vs_%d", seq_len(strata - 1L), strata)
  dimnames(Z) <- list(labels, vocab)
  list(Z = Z, labels = labels, assignment = as.character(assignment))
}
