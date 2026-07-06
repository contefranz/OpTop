if (getRversion() >= "2.15.1") {
  # docid stays declared: agg_document_stability() (deprecated, untouched)
  # calls it unqualified without importing it from quanteda
  utils::globalVariables(c("id_doc", "document", "check", "topics", ".",
                           "pval", "OpTop", "raw", "df", "pval_chisq",
                           "docid"))
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
#'   through their fitted word probabilities (see Details).
#' @param weighted_dfm A weighted `quanteda::dfm` containing word proportions
#'   for each document, built from the same counts dfm the models were
#'   fitted on; it is recommended that document ids are available via
#'   `quanteda::docid()`.
#' @param q Numeric in \eqn{(0, 1]}. Cumulative probability mass retained as
#'   "relatively important" words; the remaining words are collapsed into a
#'   single bin whose mass stays strictly below \eqn{1 - q}. Equals
#'   \eqn{1 - I^K} in the paper's notation; the default `0.95` matches the
#'   paper's numerical setting \eqn{I^K = 0.05}.
#' @param alpha Numeric in \eqn{[0, 1]}. Significance level used by the
#'   `"sequential"` and `"legacy"` selection rules (default `0.05`); ignored
#'   by `"min"`.
#' @param selection Character; how the optimal `K` is chosen (see Details):
#'   `"sequential"` (default), `"min"`, or `"legacy"` (deprecated).
#' @param calibrate Character; how the p-values are computed (see the
#'   Calibration part of Details): `"none"` (default, the asymptotic
#'   chi-square of Equation 8), `"bootstrap"` (parametric bootstrap under the
#'   fitted-model null), or `"moment"` (exact multinomial moments matched to
#'   a scaled chi-square). Not available with `selection = "legacy"`.
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
#' @param lda_models `r lifecycle::badge("deprecated")` Renamed to
#'   `topic_models`, which accepts more than `topicmodels::LDA()` fits; the
#'   old name still works but warns, and will be removed before v1.0.0.
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
#' - `"legacy"`: reproduce the pre-0.9.9 behavior exactly (rounded cutoff
#'   with the crossing word collapsed, \eqn{P_j} scaling, lower-tail p-value
#'   with 1 degree of freedom, "first `pval <= alpha`" rule). Deprecated; it
#'   only accepts `LDA_VEM` fits and will be removed before v1.0.0.
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
#' - `pval`: the p-value the selection rules use — asymptotic upper tail for
#'   `calibrate = "none"`, calibrated otherwise (legacy: the deprecated
#'   lower-tail value).
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
                          selection = c("sequential", "min", "legacy"),
                          calibrate = c("none", "bootstrap", "moment"),
                          n_boot = 200L, doc_lengths = NULL, seed = NULL,
                          n_threads = 1L, do_plot = TRUE, verbose = TRUE,
                          lda_models = lifecycle::deprecated()) {

  if (lifecycle::is_present(lda_models)) {
    lifecycle::deprecate_warn(when = "0.11.0",
                              what = "optimal_topic(lda_models)",
                              with = "optimal_topic(topic_models)")
    if (missing(topic_models)) {
      topic_models <- lda_models
    }
  }

  if (!is.list(topic_models)) {
    stop("topic_models must be a list of fitted topic models")
  }
  if (length(topic_models) == 1L) {
    stop(paste("length(topic_models) = 1.",
               "This is strange since the test should be performed",
               "on multiple topic models."))
  }
  if (!quanteda::is.dfm(weighted_dfm)) {
    stop("weighted_dfm must be a dfm")
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

  legacy <- selection == "legacy"
  if (legacy) {
    if (calibrate != "none") {
      stop('calibration is not available with selection = "legacy"')
    }
    lifecycle::deprecate_warn(
      when = "0.9.9",
      what = I('optimal_topic(selection = "legacy")'),
      details = paste("The legacy rule (lower-tail p-value with 1 degree of",
                      "freedom) predates the calibration to the published",
                      "test and will be removed before v1.0.0.")
    )
    if (!all(sapply(topic_models, is.LDA_VEM))) {
      stop(paste('selection = "legacy" reproduces the pre-0.9.9 pipeline,',
                 "which supports only LDA_VEM objects as computed by",
                 "topicmodels::LDA()"))
    }
  }

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

  if (legacy) {
    # frozen pre-0.9.9 semantics: membership is checked against the first
    # model only and one map is shared by the grid (VEM fits all see the
    # documents in input order)
    doc_check <- docs %in% topic_models[[1L]]@documents
    if (!all(doc_check)) {
      id_toremove <- which(!doc_check)
      if (length(id_toremove) < length(doc_check)) {
        cli::cli_alert_warning(
          "Removed {length(id_toremove)} document{?s} not present in the models"
        )
        weighted_dfm <- weighted_dfm[-id_toremove, ]
        docs <- docs[-id_toremove]
      } else {
        stop("Document matching went really wrong. Check docs in both weighted_dfm and in LDA@documents")
      }
    }
    doc_map <- match(docs, topic_models[[1L]]@documents) - 1L
    ks <- vapply(topic_models, function(m) as.integer(m@k), integer(1L))
  } else {
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
      tp <- optop_as_theta_phi(topic_models[[i_mod]])
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
  if (!legacy) {
    # transpose once: a document becomes a contiguous sparse column for the
    # C++ core, and every per-model call shares the same copy
    dfm_t <- Matrix::t(methods::as(weighted_dfm, "dgCMatrix"))
  }
  for (i_mod in seq_len(n_models)) {
    if (legacy) {
      core_out <- optimal_topic_core_legacy(topic_models[i_mod], weighted_dfm,
                                            q, doc_map)
    } else {
      tp <- tp_cache[[i_mod]]
      if (is.null(tp)) {
        tp <- optop_as_theta_phi(topic_models[[i_mod]])
      }
      doc_map_i <- match(docs, meta[[i_mod]]$docs) - 1L
      core_out <- optimal_topic_core(tp$theta, tp$phi, dfm_t, q, doc_map_i,
                                     calibrating, n_threads)
    }
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

  Chi_K <- data.table::as.data.table(do.call(rbind, Chi_K_rows))
  data.table::setnames(Chi_K, old = names(Chi_K), c("topic", "raw", "df"))
  Chi_K[, OpTop := raw / df]
  if (legacy) {
    # pre-0.9.9 convention: lower tail of the standardized statistic on 1 df
    Chi_K[, pval := stats::pchisq(OpTop, df = 1L, lower.tail = TRUE)]
  } else {
    # Eq. (8): the raw statistic is chi-square with df = sum_j P_j under
    # adequacy, and misspecification only inflates it, so the p-value is the
    # upper tail
    Chi_K[, pval := stats::pchisq(raw, df = df, lower.tail = FALSE)]
  }
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
  } else if (selection == "min") {
    best_topic <- global_min
    rule <- "global minimum"
  } else {
    alpha_min <- Chi_K[pval <= alpha][1L]
    if (alpha == 0 || all(is.na(alpha_min))) {
      best_topic <- global_min
      rule <- "global minimum (legacy)"
    } else if (global_min$topic > alpha_min$topic) {
      best_topic <- alpha_min
      rule <- paste0("legacy significance level of ", alpha)
    } else {
      best_topic <- global_min
      rule <- "global minimum (legacy)"
    }
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
