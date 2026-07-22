if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("id_doc", "document", "check", "topic", "topics",
                           ".", "pval", "OpTop", "raw", "df", "pval_chisq"))
}
#' Find the optimal number of topics from a pool of topic models
#'
#' Identify the number of topics that best describes the corpus with the
#' Test 1 statistic of Lewis and Grossetti (2022), a Pearson chi-square
#' goodness-of-fit test on each document's word distribution. The routine
#' evaluates each model of the grid, optionally calibrates the p-values under
#' the fitted-model null, and selects the optimal topic count with one of
#' three rules; see Details.
#'
#' @param topic_models A list of fitted topic models spanning the candidate
#'   values of \eqn{K}, one model per \eqn{K} (an unordered grid is sorted by
#'   increasing \eqn{K}). Supported classes: `topicmodels::LDA()` fits (VEM
#'   or Gibbs), `topicmodels::CTM()` fits, the seededlda models
#'   (`seededlda::textmodel_lda()`, `seededlda::textmodel_seededlda()`,
#'   `seededlda::textmodel_seqlda()`) and NLPstudio fits (`nlp_topic_fit`).
#'   Classes can be mixed within a grid as long as every model was fitted on
#'   the same corpus and vocabulary; the test consumes the models only
#'   through their fitted word probabilities (see Details). Elements may
#'   also be loader functions returning a fit: each is materialized on
#'   demand and released, so very large grids never sit in memory at once
#'   (extracted weights are still cached under a fixed budget when they
#'   fit).
#' @param weighted_dfm A weighted `quanteda::dfm` containing word proportions
#'   for each document, built from the same counts dfm the models were
#'   fitted on; it is recommended that document ids are available via
#'   `quanteda::docid()`. May also be an [optop_corpus()] of proportion
#'   shards for corpora past the size of a single dfm: shards are evaluated
#'   one at a time and combine exactly (a one-shard corpus reproduces the
#'   plain call bit for bit, calibrated p-values included; multi-shard runs
#'   agree up to floating-point summation order). With a corpus,
#'   `doc_lengths` must be named.
#' @param q Numeric in \eqn{(0, 1]}. Cumulative probability mass retained as
#'   "relatively important" words; the remaining words are collapsed into a
#'   single bin whose mass stays strictly below \eqn{1 - q}. Equals
#'   \eqn{1 - I^K} in the paper's notation; the default `0.95` matches the
#'   paper's numerical setting \eqn{I^K = 0.05}.
#' @param alpha Numeric in \eqn{[0, 1]}. Significance level used by the
#'   `"sequential"` selection rule (default `0.05`); ignored by `"min"`.
#' @param selection Character; how the optimal `K` is chosen (see Details):
#'   `"sequential"` (default) or `"min"`.
#' @param calibrate Character; how the p-values are computed (see the
#'   Calibration part of Details): `"none"` (default, the asymptotic
#'   chi-square of Equation 8), `"bootstrap"` (parametric bootstrap under the
#'   fitted-model null), or `"moment"` (exact multinomial moments matched to
#'   a scaled chi-square).
#' @param n_boot Integer; number of bootstrap replicates per model when
#'   `calibrate = "bootstrap"` (default `200L`, minimum `20`).
#' @param doc_lengths Numeric vector with the length (total token count) of
#'   each document; required when `calibrate != "none"`, typically
#'   `quanteda::ntoken()` of the counts dfm the models were fitted on. If
#'   named, names are matched to the document ids of `weighted_dfm`; if
#'   unnamed, the vector must follow the row order of `weighted_dfm`.
#' @param seed Optional integer passed to [set.seed()] before the bootstrap,
#'   for reproducible calibrated p-values (default `NULL`).
#' @param n_threads Integer; number of OpenMP threads used by the compiled
#'   cores (the per-document statistic loop and the bootstrap), default
#'   `1L`. Results are **bit-identical for any value**: every document owns a
#'   deterministic RNG stream and all reductions run in a fixed order, so
#'   `n_threads` only affects wall time. On builds without OpenMP (e.g. the
#'   default macOS toolchain) the value is ignored and the cores run
#'   single-threaded.
#' @param do_plot Logical; if `TRUE`, plot the standardized statistic versus
#'   topics with vertical and horizontal guides at the selected optimum and a
#'   subtitle reporting the selection method, the selected \eqn{K}, and the
#'   calibration if any (default `TRUE`).
#' @param verbose Logical; if `TRUE` (default), report progress (a `cli`
#'   progress bar across the model grid) and a selection summary. Regardless
#'   of `verbose`, corrective actions (dropping documents the models never
#'   saw, reordering models or features, the sequential rule's fallback) are
#'   always signalled.
#'
#' @details
#' For each document \eqn{j}, the fitted word probabilities are sorted in
#' decreasing order and the smallest head whose cumulative mass exceeds
#' \eqn{q} is kept (\eqn{P_j} words); the remaining words are collapsed into
#' a single "min" bin. The document contributes \eqn{P_j + 1} times its
#' Pearson term over the \eqn{P_j + 1} bins, and the corpus statistic is the
#' sum over documents (Equation 8 of the paper), asymptotically chi-square
#' with \eqn{\sum_j P_j}{sum_j P_j} degrees of freedom. The returned `OpTop`
#' column reports the *standardized* statistic (raw statistic divided by its
#' degrees of freedom, the version plotted in the paper's Figure 2), while
#' `pval` is the p-value the selection rules consume — the upper tail of the
#' raw statistic on its full degrees of freedom by default, or the calibrated
#' value when `calibrate != "none"`.
#'
#' The statistic touches each model only through its fitted word
#' probabilities \eqn{I_j = \theta_j^\top \Phi}{I_j = theta_j' Phi}, where
#' \eqn{\theta_j}{theta_j} is the document's topic distribution and
#' \eqn{\Phi}{Phi} the topic-word matrix. Any implementation that provides
#' those two (row-stochastic) matrices is therefore admissible, whatever its
#' estimation method — which is what the supported classes listed under
#' `topic_models` have in common.
#'
#' **Selection rules.**
#' - `"sequential"` (default): scan \eqn{K} upward and select the smallest
#'   \eqn{K} the test fails to reject (`pval > alpha`) — the classical
#'   sequential scheme for model order. If every model is rejected, the rule
#'   falls back to the global minimum with a warning.
#' - `"min"`: select the \eqn{K} with the minimum standardized statistic —
#'   the rule used in the published case study.
#'
#' The pre-0.9.9 `"legacy"` rule, deprecated since 0.9.9, was removed in
#' 0.16.0.
#'
#' **Calibration.**
#' The chi-square reference of Equation 8 is a yardstick rather than an exact
#' null law, for two reasons. First, the classical Pearson asymptotics hold
#' for the *count* statistic, whose scale factor is the document length
#' \eqn{N_j}; Equation 8 works on proportions scaled by the bin count
#' \eqn{P_j + 1} instead, so the statistic's null magnitude is off by roughly
#' \eqn{N_j / (P_j + 1)} per document. Second, the expected probabilities are
#' estimated from the same data they are tested against. The practical
#' consequence is that with \eqn{\sum_j P_j}{sum_j P_j} degrees of freedom in
#' the thousands the chi-square quantiles are razor-thin and upper-tail
#' p-values saturate at 0 or 1 unless the fit is genuinely borderline —
#' `alpha` is not a true Type-I error rate.
#'
#' Calibration replaces the chi-square reference with the distribution of the
#' statistic under the **conditional fitted-model null**: document \eqn{j} is
#' \eqn{\mathrm{Multinomial}(N_j, I_j)}{Multinomial(N_j, I_j)}, where
#' \eqn{I_j} are the \eqn{K}-model's fitted word probabilities and \eqn{N_j}
#' the observed document lengths (hence `doc_lengths`). Two properties make
#' this exact and fast:
#' - the per-document envelope (sorted fitted probabilities, cutoff at
#'   \eqn{q}, collapsed min bin) depends only on the model, never on the
#'   data, and the statistic touches the data only through sums over those
#'   fixed bins;
#' - a multinomial collapsed over bins is multinomial on the collapsed
#'   probabilities, so the null can be simulated *directly on the*
#'   \eqn{P_j + 1} *bins* — exactly equivalent to simulating whole documents
#'   over the vocabulary, at a tiny fraction of the cost.
#'
#' `calibrate = "bootstrap"` draws \eqn{B} = `n_boot` null replicates
#' \eqn{T^\ast}{T*} of the statistic \eqn{T} this way and reports the
#' empirical upper-tail p-value with the standard finite-sample correction,
#' \eqn{(1 + \#\{T^\ast \ge T\}) / (B + 1)}{(1 + #{T* >= T}) / (B + 1)}.
#' This is the reference method: `alpha` becomes a genuine Type-I error rate
#' with respect to the conditional null, at bootstrap resolution
#' \eqn{1 / (B + 1)}.
#'
#' `calibrate = "moment"` uses the exact multinomial moments of the
#' per-document Pearson term (Haldane, 1937): over \eqn{k} bins with
#' probabilities \eqn{p_b} and length \eqn{n}, the count statistic
#' \eqn{X^2} has
#' \eqn{E[X^2] = k - 1}{E[X^2] = k - 1} and
#' \eqn{\mathrm{Var}[X^2] = 2(k - 1) + (\sum_b 1/p_b - k^2 - 2k + 2)/n}{Var[X^2] = 2(k - 1) + (sum_b 1/p_b - k^2 - 2k + 2)/n};
#' scaling by the statistic's \eqn{k/n} factor and summing over independent
#' documents gives the null mean \eqn{\mu}{mu} and variance
#' \eqn{\sigma^2}{sigma^2} of the corpus statistic, which are matched to a
#' scaled chi-square \eqn{a\,\chi^2_\nu}{a * chisq(nu)} (Satterthwaite:
#' \eqn{a = \sigma^2 / (2\mu)}{a = sigma^2 / (2 mu)},
#' \eqn{\nu = 2\mu^2 / \sigma^2}{nu = 2 mu^2 / sigma^2}). Closed form, no
#' simulation — a fast approximation that corrects the location and scale of
#' the reference but not its higher moments.
#'
#' One caveat applies to both: the null holds the estimated
#' \eqn{\theta}{theta} and \eqn{\phi}{phi} fixed (no per-replicate re-fitting
#' of the model — the "double bootstrap" would be exact but computationally
#' prohibitive). Calibrated p-values are therefore conditional on the fitted
#' parameters and do not account for estimation noise in \eqn{\theta}{theta}
#' and \eqn{\phi}{phi}.
#'
#' **Input alignment.** `weighted_dfm` must be a `quanteda::dfm` of word
#' proportions (row-wise), built from the same counts dfm the models were
#' fitted on. Its features must match the models' common vocabulary: when
#' they are the same set in a different order the columns are reordered
#' automatically (signalled), any other mismatch is an error. Document
#' identifiers are taken from `quanteda::docid(weighted_dfm)` and matched
#' against each model's own document identifiers: documents not present in
#' *every* model are dropped (with a warning), and each retained dfm row is
#' paired with the corresponding row of \eqn{\theta}{theta} *per model*.
#' Neither the row order of `weighted_dfm` nor the order in which each model
#' saw the documents matters — alignment is always by identifier.
#'
#' **Performance note.** Both hot paths run in C++ compiled code and
#' parallelize over documents with OpenMP (see `n_threads`): the statistic
#' evaluates each block of documents with one BLAS product and sorts the
#' fitted probabilities per document, and the bootstrap draws the null
#' replicates directly on the collapsed bins, fused with the Pearson
#' reduction, at about \eqn{B \times \mathrm{df}}{B x df} floating-point
#' operations per model. See the "Computational efficiency" section of the
#' vignette for measurements.
#'
#' @return A `data.table` with columns:
#' - `topic`: integer number of topics (\eqn{K}).
#' - `OpTop`: standardized Test 1 statistic (raw statistic divided by `df`).
#' - `df`: degrees of freedom \eqn{\sum_j P_j}{sum_j P_j} of the raw
#'   statistic.
#' - `pval`: the p-value the selection rules use: asymptotic upper tail for
#'   `calibrate = "none"`, calibrated otherwise.
#' - `pval_chisq`: only when `calibrate != "none"` — the uncalibrated
#'   asymptotic chi-square p-value, for comparison.
#'
#' @examples
#' \dontrun{
#' # Asymptotic p-values (Equation 8)
#' test1 <- optimal_topic(topic_models = lda_list,
#'                        weighted_dfm = weighted_dfm,
#'                        q = 0.95,
#'                        alpha = 0.05,
#'                        selection = "sequential")
#'
#' # Bootstrap-calibrated p-values: document lengths come from the counts dfm
#' test1_cal <- optimal_topic(topic_models = lda_list,
#'                            weighted_dfm = weighted_dfm,
#'                            calibrate = "bootstrap",
#'                            n_boot = 200,
#'                            doc_lengths = quanteda::ntoken(counts_dfm),
#'                            seed = 42)
#' }
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2022). A statistical approach for optimal
#' topic model identification. *Journal of Machine Learning Research*,
#' 23(58), 1--20. <https://jmlr.org/papers/v23/19-297.html>
#'
#' Haldane, J. B. S. (1937). The exact value of the moments of the
#' distribution of chi-square. *Biometrika*, 29, 133--143.
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin*, 2(6), 110--114.
#'
#' Davison, A. C. and Hinkley, D. V. (1997). *Bootstrap Methods and their
#' Application*. Cambridge University Press, Cambridge.
#'
#' @seealso [topicmodels::LDA()], [topicmodels::CTM()],
#'   [seededlda::textmodel_lda()]
#'
#' @import data.table
#' @export

optimal_topic <- function(topic_models, weighted_dfm, q = 0.95, alpha = 0.05,
                          selection = c("sequential", "min"),
                          calibrate = c("none", "bootstrap", "moment"),
                          n_boot = 200L, doc_lengths = NULL, seed = NULL,
                          n_threads = 1L, do_plot = TRUE, verbose = TRUE) {

  if (!is.list(topic_models)) {
    stop("topic_models must be a list of fitted topic models")
  }
  if (length(topic_models) == 1L) {
    stop(paste("length(topic_models) = 1.",
               "This is strange since the test should be performed",
               "on multiple topic models."))
  }
  if (!quanteda::is.dfm(weighted_dfm) && !.optop_is_corpus(weighted_dfm)) {
    stop("weighted_dfm must be a dfm or an optop_corpus")
  }
  if (!is.numeric(q)) {
    stop("q must be a numeric")
  }
  if (!is.numeric(alpha)) {
    stop("alpha must be a numeric")
  }
  selection <- match.arg(selection)
  calibrate <- match.arg(calibrate)
  if (!is.logical(verbose)) {
    stop("verbose must be either TRUE or FALSE")
  }
  if (!is.numeric(n_threads) || length(n_threads) != 1L ||
      !is.finite(n_threads) || n_threads < 1) {
    stop("n_threads must be a single integer >= 1")
  }
  n_threads <- as.integer(n_threads)

  calibrating <- calibrate != "none"
  if (calibrating) {
    if (is.null(doc_lengths)) {
      stop(paste("calibrate != \"none\" needs doc_lengths: the total token",
                 "count of each document, typically quanteda::ntoken() of",
                 "the counts dfm the models were fitted on"))
    }
    if (!is.numeric(doc_lengths) || any(!is.finite(doc_lengths)) ||
        any(doc_lengths <= 0)) {
      stop("doc_lengths must be a numeric vector of positive document lengths")
    }
    if (!is.numeric(n_boot) || length(n_boot) != 1L || n_boot < 20) {
      stop("n_boot must be a single number >= 20")
    }
    n_boot <- as.integer(n_boot)
  }

  tic <- proc.time()
  if (verbose) {
    cli::cli_h2("Optimal topic selection")
  }
  if (.optop_is_corpus(weighted_dfm)) {
    return(.optop_ot_corpus(weighted_dfm, topic_models, q, alpha, selection,
                            calibrate, n_boot, doc_lengths, seed, n_threads,
                            do_plot, verbose, tic))
  }
  docs <- as.character(quanteda::docid(weighted_dfm))
  if (anyDuplicated(docs)) {
    stop(paste("weighted_dfm has duplicated document ids; document",
               "alignment would be ambiguous"))
  }
  if (calibrating && is.null(names(doc_lengths)) &&
      length(doc_lengths) != length(docs)) {
    stop(paste("unnamed doc_lengths must have one entry per document of",
               "weighted_dfm, in the same row order"))
  }

  # adapt every model to the common (theta, phi, K, docs, terms) contract.
  # The extracted weight matrices are kept for the evaluation loop as long
  # as the whole grid fits a fixed memory budget (typical grids amount to a
  # few MB); past the budget the remaining models are re-extracted one at a
  # time inside the loop, so peak memory stays bounded on very large
  # vocabularies
  meta <- vector("list", length(topic_models))
  tp_cache <- vector("list", length(topic_models))
  cache_bytes <- 0
  cache_budget <- 256 * 1024^2
  for (i_mod in seq_along(topic_models)) {
    tp <- optop_as_theta_phi(.optop_materialize_model(topic_models[[i_mod]]))
    .optop_validate_theta_phi(tp, class(topic_models[[i_mod]]))
    meta[[i_mod]] <- tp[c("K", "docs", "terms")]
    cache_bytes <- cache_bytes + 8 * (length(tp$theta) + length(tp$phi))
    if (cache_bytes <= cache_budget) {
      tp_cache[[i_mod]] <- tp[c("theta", "phi")]
    }
  }
  ks <- vapply(meta, function(m) as.integer(m$K), integer(1L))
  if (anyDuplicated(ks)) {
    stop(paste("topic_models contains more than one model with the same",
               "number of topics; the grid must have one model per K"))
  }
  if (is.unsorted(ks)) {
    ord <- order(ks)
    topic_models <- topic_models[ord]
    meta <- meta[ord]
    tp_cache <- tp_cache[ord]
    ks <- ks[ord]
    cli::cli_alert_warning(
      "Reordered topic_models by increasing number of topics"
    )
  }

  # one common vocabulary across the grid, in one common order
  ref_terms <- meta[[1L]]$terms
  same_terms <- vapply(meta[-1L], function(m) identical(m$terms, ref_terms),
                       logical(1L))
  if (!all(same_terms)) {
    stop(paste("all models must share the same vocabulary in the same",
               "order; refit the grid on a common dfm"))
  }
  feats <- quanteda::featnames(weighted_dfm)
  if (!identical(feats, ref_terms)) {
    if (length(feats) == length(ref_terms) && !anyDuplicated(feats) &&
        setequal(feats, ref_terms)) {
      # a permutation preserves the row proportions, so it is safe to fix
      weighted_dfm <- weighted_dfm[, ref_terms]
      cli::cli_alert_warning(
        "Reordered the features of weighted_dfm to the models' vocabulary"
      )
    } else {
      stop(paste("the features of weighted_dfm do not match the models'",
                 "vocabulary. Rebuild the weighted dfm from the counts dfm",
                 "the models were fitted on, e.g.",
                 'quanteda::dfm_weight(counts_dfm, scheme = "prop")'))
    }
  }
  row_mass <- Matrix::rowSums(weighted_dfm)
  if (any(abs(row_mass - 1) > 1e-6)) {
    cli::cli_alert_warning(paste(
      "weighted_dfm rows do not sum to 1: expected word proportions as",
      'produced by quanteda::dfm_weight(scheme = "prop")'
    ))
  }

  # keep only the documents every model has seen: identifiers are matched
  # per model, so the whole grid is evaluated on one fixed document set and
  # the statistics stay comparable across K (and across engines)
  keep <- rep(TRUE, length(docs))
  for (m in meta) {
    keep <- keep & docs %in% m$docs
  }
  if (!all(keep)) {
    id_toremove <- which(!keep)
    if (length(id_toremove) == length(docs)) {
      stop(paste("Document matching went really wrong. Check the document",
                 "ids of weighted_dfm against those stored in the models"))
    }
    cli::cli_alert_warning(
      "Removed {length(id_toremove)} document{?s} not present in the models"
    )
    weighted_dfm <- weighted_dfm[-id_toremove, ]
    docs <- docs[-id_toremove]
    if (calibrating && is.null(names(doc_lengths))) {
      doc_lengths <- doc_lengths[-id_toremove]
    }
  }
  if (calibrating && !is.null(names(doc_lengths))) {
    doc_lengths <- doc_lengths[docs]
    if (anyNA(doc_lengths)) {
      stop("doc_lengths is missing entries for some documents of weighted_dfm")
    }
  }

  n_models <- length(topic_models)
  if (verbose) {
    cli::cli_alert_info(paste(
      "Evaluating {n_models} models on {length(docs)} document{?s} and",
      "{quanteda::nfeat(weighted_dfm)} features (q = {q}, alpha = {alpha},",
      "selection = {selection})"
    ))
    if (calibrate == "bootstrap") {
      seed_note <- if (is.null(seed)) "" else paste0(", seed = ", seed)
      cli::cli_alert_info(paste0(
        "Calibrating p-values by parametric bootstrap (B = {n_boot}",
        seed_note, ") under the fitted-model null"
      ))
    } else if (calibrate == "moment") {
      cli::cli_alert_info(paste(
        "Calibrating p-values by moment matching (Satterthwaite scaled",
        "chi-square) under the fitted-model null"
      ))
    }
    cli::cli_progress_bar("Processing model grid", total = n_models)
  }
  if (calibrating && !is.null(seed)) {
    set.seed(seed)
  }

  Chi_K_rows <- vector("list", n_models)
  pval_cal <- rep(NA_real_, n_models)
  # transpose once: a document becomes a contiguous sparse column for the
  # C++ core. The raw CSC slots (@p, @i, @x) cross the boundary as
  # zero-copy views, so no per-model conversion or copy of the corpus is
  # ever made and the J x W shape is unbounded
  dfm_t <- Matrix::t(methods::as(weighted_dfm, "dgCMatrix"))
  dfm_terms <- nrow(dfm_t)
  for (i_mod in seq_len(n_models)) {
    tp <- tp_cache[[i_mod]]
    if (is.null(tp)) {
      tp <- optop_as_theta_phi(.optop_materialize_model(topic_models[[i_mod]]))
    }
    doc_map_i <- match(docs, meta[[i_mod]]$docs) - 1L
    core_out <- optimal_topic_core(tp$theta, tp$phi, dfm_t@p, dfm_t@i,
                                   dfm_t@x, dfm_terms, q, doc_map_i,
                                   calibrating, n_threads)
    Chi_K_rows[[i_mod]] <- core_out$stat
    if (calibrating) {
      T_obs <- core_out$stat[1L, 2L]
      if (calibrate == "bootstrap") {
        # the flattened envelope goes to the compiled core as exported: no
        # per-document list, no unlist round trip
        T_null <- .optop_boot_null(doc_lengths = doc_lengths, n_boot = n_boot,
                                   n_threads = n_threads,
                                   bin_probs = core_out$bin_probs,
                                   bin_counts = core_out$bin_counts)
        pval_cal[i_mod] <- (1 + sum(T_null >= T_obs)) / (n_boot + 1)
      } else {
        probs <- .optop_split_envelope(core_out$bin_probs, core_out$bin_counts)
        mm <- .optop_moment_null(probs, doc_lengths)
        pval_cal[i_mod] <- stats::pchisq(T_obs / mm$a, df = mm$nu,
                                         lower.tail = FALSE)
      }
    }
    if (verbose) {
      cli::cli_progress_update(status = paste0("k = ", ks[i_mod]))
    }
  }
  if (verbose) {
    cli::cli_progress_done()
  }

  .optop_ot_finish(Chi_K_rows, pval_cal, calibrating,
                   calibrate, n_boot, selection, alpha, do_plot,
                   verbose, tic)
}

#' Shared assembly, selection, and reporting tail of optimal_topic()
#'
#' Consumes the per-model statistic rows (K, raw, df) plus the calibrated
#' p-values, builds the returned table, applies the selection rule, and
#' handles the plot and the verbose reporting. One implementation serves
#' the single-matrix and the sharded-corpus evaluation paths.
#'
#' @keywords internal
.optop_ot_finish <- function(Chi_K_rows, pval_cal, calibrating,
                             calibrate, n_boot, selection, alpha, do_plot,
                             verbose, tic) {
  Chi_K <- data.table::as.data.table(do.call(rbind, Chi_K_rows))
  data.table::setnames(Chi_K, old = names(Chi_K), c("topic", "raw", "df"))
  Chi_K[, OpTop := raw / df]
  # Eq. (8): the raw statistic is chi-square with df = sum_j P_j under
  # adequacy, and misspecification only inflates it, so the p-value is the
  # upper tail
  Chi_K[, pval := stats::pchisq(raw, df = df, lower.tail = FALSE)]
  Chi_K[, raw := NULL]
  if (calibrating) {
    # the selection rules read pval: the calibrated value takes its place and
    # the asymptotic one stays alongside for comparison
    Chi_K[, pval_chisq := pval]
    Chi_K[, pval := pval_cal]
    data.table::setcolorder(Chi_K, c("topic", "OpTop", "df", "pval",
                                     "pval_chisq"))
  } else {
    data.table::setcolorder(Chi_K, c("topic", "OpTop", "df", "pval"))
  }

  cal_suffix <- switch(calibrate,
                       none = "",
                       bootstrap = ", bootstrap-calibrated",
                       moment = ", moment-calibrated")

  global_min <- Chi_K[, .SD[which.min(OpTop)]]
  method_label <- selection
  if (selection == "sequential") {
    accepted <- Chi_K[pval > alpha][1L]
    if (all(is.na(accepted))) {
      cli::cli_alert_warning(paste(
        "No model passes the adequacy test at alpha = {alpha};",
        "falling back to the global minimum"
      ))
      best_topic <- global_min
      rule <- "global minimum (sequential fallback)"
      method_label <- "sequential (fallback: min)"
    } else {
      best_topic <- accepted
      rule <- paste0("sequential adequacy scan at alpha = ", alpha)
    }
  } else {
    best_topic <- global_min
    rule <- "global minimum"
  }
  if (verbose) {
    cli::cli_alert_success(
      "Optimal model has {best_topic$topic} topics (selected by {rule}{cal_suffix})"
    )
  }

  if (do_plot) {
    x_min <- best_topic$topic
    y_min <- best_topic$OpTop
    subtitle_expr <- if (calibrating) {
      cal_label <- if (calibrate == "bootstrap") {
        paste0("bootstrap (B = ", n_boot, ")")
      } else {
        "moment"
      }
      bquote("Method:" ~ .(method_label) ~ "|" ~ K[opt] == .(x_min) ~ "|" ~
               "calibrated:" ~ .(cal_label))
    } else {
      bquote("Method:" ~ .(method_label) ~ "|" ~ K[opt] == .(x_min))
    }
    p1 <- ggplot2::ggplot(Chi_K) +
      ggplot2::geom_line(ggplot2::aes(x = topic, y = OpTop),
                         linewidth = 0.8, color = "royalblue") +
      ggplot2::geom_hline(yintercept = y_min, color = "black", linetype = 2L) +
      ggplot2::geom_vline(xintercept = x_min, color = "black", linetype = 2L) +
      ggplot2::annotate("point", x = x_min, y = y_min,
                        color = "red", shape = 4L, size = 4L) +
      ggplot2::xlab("Topics") + ggplot2::ylab(expression(OpTop[J]^{"K"})) +
      ggplot2::ggtitle("Optimal Topic Plot") +
      ggplot2::labs(subtitle = subtitle_expr) +
      ggplot2::theme_bw()
    print(p1)
  }

  toc <- proc.time()
  runtime <- toc - tic
  if (verbose) {
    cli::cli_alert_info("Completed in {round(runtime[3L], 2)}s")
  }
  return(Chi_K)
}

#' Sharded evaluation path of optimal_topic()
#'
#' Evaluates the model grid over an [optop_corpus()], one shard at a time:
#' per shard and per model, the compiled core contributes the shard's raw
#' statistic and degrees of freedom, which add exactly across shards
#' because the Test 1 statistic is a sum of independent per-document terms.
#' Calibration composes the same way: bootstrap replicates are drawn per
#' shard with the document streams keyed by the global document index
#' (`doc_offset`), and the closed-form moments accumulate their per-shard
#' mean and variance contributions. Models supplied as loader functions are
#' materialized on demand; extracted weight matrices are cached across
#' shards under the usual fixed memory budget, so cheap grids avoid
#' reloading while very large grids stream one model at a time.
#'
#' One requirement is specific to the sharded path: with calibration,
#' `doc_lengths` must be named (row order across shards is not a reliable
#' identity at this scale).
#'
#' @keywords internal
.optop_ot_corpus <- function(corpus, topic_models, q, alpha, selection,
                             calibrate, n_boot, doc_lengths, seed, n_threads,
                             do_plot, verbose, tic) {
  calibrating <- calibrate != "none"
  if (calibrating && is.null(names(doc_lengths))) {
    stop(paste("corpus input needs named doc_lengths: names are matched to",
               "the document ids of every shard"))
  }

  # meta pass: materialize each model once for K, vocabulary, and document
  # ids; identical document vectors are stored once (the common case)
  n_models <- length(topic_models)
  metas <- vector("list", n_models)
  docs_ref <- NULL
  for (i in seq_len(n_models)) {
    tp <- optop_as_theta_phi(.optop_materialize_model(topic_models[[i]]))
    .optop_validate_theta_phi(tp, class(topic_models[[i]]))
    own_docs <- if (is.null(docs_ref)) {
      docs_ref <- tp$docs
      NULL
    } else if (identical(tp$docs, docs_ref)) {
      NULL
    } else {
      tp$docs
    }
    metas[[i]] <- list(K = tp$K, terms = tp$terms, docs = own_docs)
  }
  ks <- vapply(metas, function(m) as.integer(m$K), integer(1L))
  if (anyDuplicated(ks)) {
    stop(paste("topic_models contains more than one model with the same",
               "number of topics; the grid must have one model per K"))
  }
  if (is.unsorted(ks)) {
    ord <- order(ks)
    topic_models <- topic_models[ord]
    metas <- metas[ord]
    ks <- ks[ord]
    cli::cli_alert_warning(
      "Reordered topic_models by increasing number of topics"
    )
  }
  ref_terms <- metas[[1L]]$terms
  same_terms <- vapply(metas[-1L], function(m) identical(m$terms, ref_terms),
                       logical(1L))
  if (!all(same_terms)) {
    stop(paste("all models must share the same vocabulary in the same",
               "order; refit the grid on a common dfm"))
  }

  # extracted weight matrices are cached across shards under the fixed
  # budget; past it, models are re-materialized per shard
  cache_budget <- 256 * 1024^2
  tp_cache <- vector("list", n_models)
  cache_bytes <- 0
  tp_get <- function(i) {
    if (!is.null(tp_cache[[i]])) {
      return(tp_cache[[i]])
    }
    tp <- optop_as_theta_phi(.optop_materialize_model(topic_models[[i]]))
    bytes <- 8 * (length(tp$theta) + length(tp$phi))
    if (cache_bytes + bytes <= cache_budget) {
      tp_cache[[i]] <<- tp[c("theta", "phi", "docs")]
      cache_bytes <<- cache_bytes + bytes
      return(tp_cache[[i]])
    }
    tp[c("theta", "phi", "docs")]
  }

  n_sh <- corpus$n_shards
  if (verbose) {
    cli::cli_alert_info(paste(
      "Evaluating {n_models} models over a corpus of {n_sh} shard{?s}",
      "({length(ref_terms)} features; q = {q}, alpha = {alpha},",
      "selection = {selection})"
    ))
    if (calibrate == "bootstrap") {
      seed_note <- if (is.null(seed)) "" else paste0(", seed = ", seed)
      cli::cli_alert_info(paste0(
        "Calibrating p-values by parametric bootstrap (B = {n_boot}",
        seed_note, ") under the fitted-model null"
      ))
    } else if (calibrate == "moment") {
      cli::cli_alert_info(paste(
        "Calibrating p-values by moment matching (Satterthwaite scaled",
        "chi-square) under the fitted-model null"
      ))
    }
    cli::cli_progress_bar("Processing corpus shards",
                          total = n_sh * n_models)
  }

  # one bootstrap seed per model, drawn upfront so every shard of a model
  # shares the stream base (set.seed() governs reproducibility as usual)
  seeds_model <- NULL
  if (calibrate == "bootstrap") {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    seeds_model <- sample.int(.Machine$integer.max, n_models)
  }

  raw_sum <- numeric(n_models)
  df_sum <- numeric(n_models)
  T_null_sum <- if (calibrate == "bootstrap") {
    matrix(0, n_boot, n_models)
  } else {
    NULL
  }
  mm_mu <- numeric(n_models)
  mm_s2 <- numeric(n_models)
  T_kept <- 0
  n_dropped <- 0
  warned_perm <- FALSE
  warned_mass <- FALSE

  for (s in seq_len(n_sh)) {
    m <- .optop_corpus_shard(corpus, s)
    docs_s <- rownames(m)
    if (is.null(docs_s) || anyDuplicated(docs_s)) {
      stop(sprintf("shard %d needs unique document identifiers", s))
    }
    feats <- colnames(m)
    if (!identical(feats, ref_terms)) {
      if (length(feats) == length(ref_terms) && !anyDuplicated(feats) &&
          setequal(feats, ref_terms)) {
        m <- m[, ref_terms]
        if (!warned_perm) {
          cli::cli_alert_warning(
            "Reordered shard features to the models' vocabulary"
          )
          warned_perm <- TRUE
        }
      } else {
        stop(sprintf(paste("the features of shard %d do not match the",
                           "models' vocabulary"), s))
      }
    }
    if (!warned_mass) {
      row_mass <- Matrix::rowSums(m)
      if (any(abs(row_mass - 1) > 1e-6)) {
        cli::cli_alert_warning(paste(
          "corpus rows do not sum to 1: expected word proportions as",
          'produced by quanteda::dfm_weight(scheme = "prop")'
        ))
        warned_mass <- TRUE
      }
    }

    keep <- rep(TRUE, length(docs_s))
    keep <- keep & docs_s %in% docs_ref
    for (mt in metas) {
      if (!is.null(mt$docs)) {
        keep <- keep & docs_s %in% mt$docs
      }
    }
    if (!all(keep)) {
      n_dropped <- n_dropped + sum(!keep)
      m <- m[keep, , drop = FALSE]
      docs_s <- docs_s[keep]
    }
    if (!length(docs_s)) {
      for (i in seq_len(n_models)) {
        if (verbose) cli::cli_progress_update()
      }
      next
    }

    dl_s <- NULL
    if (calibrating) {
      dl_s <- doc_lengths[docs_s]
      if (anyNA(dl_s)) {
        stop("doc_lengths is missing entries for some corpus documents")
      }
    }

    dfm_t <- Matrix::t(m)
    dfm_terms <- nrow(dfm_t)

    for (i in seq_len(n_models)) {
      tp <- tp_get(i)
      model_docs <- if (is.null(metas[[i]]$docs)) docs_ref else
        metas[[i]]$docs
      doc_map_i <- match(docs_s, model_docs) - 1L
      core_out <- optimal_topic_core(tp$theta, tp$phi, dfm_t@p, dfm_t@i,
                                     dfm_t@x, dfm_terms, q, doc_map_i,
                                     calibrating, n_threads)
      raw_sum[i] <- raw_sum[i] + core_out$stat[1L, 2L]
      df_sum[i] <- df_sum[i] + core_out$stat[1L, 3L]
      if (calibrate == "bootstrap") {
        T_null_sum[, i] <- T_null_sum[, i] +
          .optop_boot_null(doc_lengths = dl_s, n_boot = n_boot,
                           seed = seeds_model[i], n_threads = n_threads,
                           bin_probs = core_out$bin_probs,
                           bin_counts = core_out$bin_counts,
                           doc_offset = T_kept)
      } else if (calibrate == "moment") {
        probs <- .optop_split_envelope(core_out$bin_probs,
                                       core_out$bin_counts)
        mm <- .optop_moment_null(probs, dl_s)
        mm_mu[i] <- mm_mu[i] + mm$mu
        mm_s2[i] <- mm_s2[i] + mm$sigma2
      }
      if (verbose) {
        cli::cli_progress_update(
          status = paste0("shard ", s, "/", n_sh, ", k = ", ks[i])
        )
      }
    }
    T_kept <- T_kept + length(docs_s)
  }
  if (verbose) {
    cli::cli_progress_done()
  }
  if (n_dropped > 0) {
    cli::cli_alert_warning(
      "Removed {n_dropped} document{?s} not present in the models"
    )
  }
  if (T_kept == 0) {
    stop(paste("Document matching went really wrong. Check the document",
               "ids of the corpus against those stored in the models"))
  }
  if (verbose) {
    cli::cli_alert_info(
      "Evaluated {T_kept} document{?s} across {n_sh} shard{?s}"
    )
  }

  pval_cal <- rep(NA_real_, n_models)
  if (calibrate == "bootstrap") {
    pval_cal <- vapply(seq_len(n_models), function(i) {
      (1 + sum(T_null_sum[, i] >= raw_sum[i])) / (n_boot + 1)
    }, numeric(1))
  } else if (calibrate == "moment") {
    pval_cal <- vapply(seq_len(n_models), function(i) {
      a <- mm_s2[i] / (2 * mm_mu[i])
      nu <- 2 * mm_mu[i]^2 / mm_s2[i]
      stats::pchisq(raw_sum[i] / a, df = nu, lower.tail = FALSE)
    }, numeric(1))
  }

  Chi_K_rows <- lapply(seq_len(n_models), function(i) {
    matrix(c(ks[i], raw_sum[i], df_sum[i]), 1L, 3L)
  })
  .optop_ot_finish(Chi_K_rows, pval_cal,
                   calibrating = calibrating, calibrate = calibrate,
                   n_boot = n_boot, selection = selection, alpha = alpha,
                   do_plot = do_plot, verbose = verbose, tic = tic)
}
