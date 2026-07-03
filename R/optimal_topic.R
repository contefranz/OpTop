if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("id_word", "id_doc", "weighted_dfm",
                           "word_prop", "document",
                           "word_sum", "check", "topics",
                           ".", "chisquare", "chisquare_mod",
                           "row_cut", "chi_sum", "word_prop_hat",
                           "word_prop_hat_cum", "pval", "docid",
                           "OpTop", "raw", "df"))
}
#' Find the optimal number of topics from a pool of LDA models
#'
#' Identify the number of topics that best describes the corpus with the
#' Test 1 statistic of Lewis and Grossetti (2022), a Pearson chi-square
#' goodness-of-fit test on each document's word distribution. The routine
#' evaluates each model of the grid and selects the optimal topic count with
#' one of three rules; see Details.
#'
#' @param lda_models A list of `topicmodels::LDA` objects (VEM), ordered by
#'   increasing number of topics. The grid should span the candidate values of `K`.
#' @param weighted_dfm A weighted `quanteda::dfm` containing word proportions
#'   for each document; it is recommended that document ids are available via
#'   `quanteda::docid()`.
#' @param q Numeric in `(0, 1]`. Cumulative probability mass retained as
#'   "relatively important" words; the remaining words are collapsed into a
#'   single bin whose mass stays strictly below `1 - q`. Equals `1 - I^K` in
#'   the paper's notation; the default `0.95` matches the paper's numerical
#'   setting `I^K = 0.05`.
#' @param alpha Numeric in `[0, 1]`. Significance level used by the
#'   `"sequential"` and `"legacy"` selection rules (default `0.05`); ignored
#'   by `"min"`.
#' @param selection Character; how the optimal `K` is chosen (see Details):
#'   `"sequential"` (default), `"min"`, or `"legacy"` (deprecated).
#' @param do_plot Logical; if `TRUE`, plot the standardized statistic versus
#'   topics with vertical and horizontal guides at the selected optimum
#'   (default `TRUE`).
#' @param verbose Logical; if `TRUE`, report progress (a `cli` progress bar
#'   across the model grid) and a selection summary (default `FALSE`).
#'   Regardless of `verbose`, dropping documents that the models never saw
#'   and the sequential rule's fallback are always signalled.
#'
#' @details
#' For each document, the fitted word probabilities are sorted in decreasing
#' order and the smallest head whose cumulative mass exceeds `q` is kept
#' (`P_j` words); the remaining words are collapsed into a single "min" bin.
#' The document contributes `(P_j + 1)` times its Pearson term over the
#' `P_j + 1` bins, and the corpus statistic is the sum over documents
#' (Equation 8 of the paper), asymptotically chi-square with
#' `df = sum_j P_j` degrees of freedom. The returned `OpTop` column reports
#' the *standardized* statistic (raw statistic divided by `df`, the version
#' plotted in the paper's Figure 2), while `pval` is the upper-tail p-value
#' of the raw statistic on its full degrees of freedom.
#'
#' **Selection rules.**
#' - `"sequential"` (default): scan `K` upward and select the smallest `K`
#'   the test fails to reject (`pval > alpha`) — the classical sequential
#'   scheme for model order. If every model is rejected, the rule falls back
#'   to the global minimum with a warning.
#' - `"min"`: select the `K` with the minimum standardized statistic — the
#'   rule used in the published case study.
#' - `"legacy"`: reproduce the pre-0.9.9 behavior exactly (rounded cutoff
#'   with the crossing word collapsed, `P_j` scaling, lower-tail p-value with
#'   1 degree of freedom, "first `pval <= alpha`" rule). Deprecated; it will
#'   be removed before v1.0.0.
#'
#' **A calibration caveat.** With `df = sum_j P_j` in the thousands, the
#' chi-square reference distribution is extremely concentrated, so upper-tail
#' p-values tend to saturate at 0 or 1 unless the fit is genuinely borderline.
#' When the sequential rule degenerates (all 0 or all 1), the shape of the
#' standardized statistic — and its minimum — carries the information, which
#' is why the paper's case study uses the `"min"` rule.
#'
#' **Input alignment.** `weighted_dfm` must be a `quanteda::dfm` of word
#' proportions (row-wise). Document identifiers are taken from
#' `quanteda::docid(weighted_dfm)` and matched to the LDA fits; documents not
#' present in the first model’s `@documents` slot are dropped (with a warning)
#' to ensure alignment. The remaining documents are paired with the rows of
#' `@gamma` by identifier, so the row order of `weighted_dfm` does not need to
#' follow the order in which the models saw the documents.
#'
#' **Performance note.** The core computation is delegated to C++ compiled code
#' to handle high-dimensional vocabularies efficiently.
#'
#' @return A `data.table` with columns:
#' - `topic`: integer number of topics (`K`).
#' - `OpTop`: standardized Test 1 statistic (raw statistic / `df`).
#' - `df`: degrees of freedom `sum_j P_j` of the raw statistic.
#' - `pval`: p-value associated with the raw statistic (upper tail for
#'   `"sequential"`/`"min"`; the legacy lower-tail value for `"legacy"`).
#'
#' @examples
#' \dontrun{
#' # Compute word proportions from a corpus objects
#' test1 <- optimal_topic(lda_models = lda_list,
#'                        weighted_dfm = weighted_dfm,
#'                        q = 0.95,
#'                        alpha = 0.05,
#'                        selection = "sequential",
#'                        verbose = TRUE)
#' }
#'
#' @references
#' Lewis, C. M. and Grossetti, F. (2022). A statistical approach for optimal
#' topic model identification. *Journal of Machine Learning Research*,
#' 23(58), 1--20. <https://jmlr.org/papers/v23/19-297.html>
#'
#' @seealso [topicmodels::LDA()]
#'
#' @import data.table
#' @export

optimal_topic <- function(lda_models, weighted_dfm, q = 0.95, alpha = 0.05,
                          selection = c("sequential", "min", "legacy"),
                          do_plot = TRUE, verbose = FALSE) {

  if (!is.list(lda_models)) {
    stop("lda_models must be a list")
  }
  if (length(lda_models) == 1L) {
    stop(paste("length(lda_models) = 1.",
               "This is strange since the test should be perfomed",
               "on multiple LDA models."))
  }
  if (!all(sapply(lda_models, is.LDA_VEM))) {
    stop(paste("lda_models must contain LDA_VEM obects as computed",
               "by topicmodels::LDA()"))
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
  if (!is.logical(verbose)) {
    stop("verbose must be either TRUE or FALSE")
  }

  if (selection == "legacy") {
    lifecycle::deprecate_warn(
      when = "0.9.9",
      what = I('optimal_topic(selection = "legacy")'),
      details = paste("The legacy rule (lower-tail p-value with 1 degree of",
                      "freedom) predates the calibration to the published",
                      "test and will be removed before v1.0.0.")
    )
  }

  tic <- proc.time()
  if (verbose) {
    cli::cli_h2("Optimal topic selection")
  }
  docs <- as.character(quanteda::docid(weighted_dfm))

  # drop documents the models never saw; all models are assumed to share the
  # document set of the first one (a document the LDA drops is dropped for
  # every k), so the check runs against lda_models[[1L]] only
  doc_check <- docs %in% lda_models[[1L]]@documents
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

  # map each dfm row to the corresponding row of @gamma (0-based for C++);
  # membership was checked above, so no NA can survive the match
  doc_map <- match(docs, lda_models[[1L]]@documents) - 1L

  n_models <- length(lda_models)
  if (verbose) {
    cli::cli_alert_info(paste(
      "Evaluating {n_models} models on {length(docs)} document{?s} and",
      "{quanteda::nfeat(weighted_dfm)} features (q = {q}, alpha = {alpha},",
      "selection = {selection})"
    ))
    cli::cli_progress_bar("Processing LDA grid", total = n_models)
  }

  legacy <- selection == "legacy"
  Chi_K_rows <- vector("list", n_models)
  for (i_mod in seq_len(n_models)) {
    Chi_K_rows[[i_mod]] <- optimal_topic_core(lda_models[i_mod],
                                              weighted_dfm, q, doc_map,
                                              legacy)
    if (verbose) {
      cli::cli_progress_update(
        status = paste0("k = ", lda_models[[i_mod]]@k)
      )
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
  data.table::setcolorder(Chi_K, c("topic", "OpTop", "df", "pval"))

  global_min <- Chi_K[, .SD[which.min(OpTop)]]
  if (selection == "sequential") {
    accepted <- Chi_K[pval > alpha][1L]
    if (all(is.na(accepted))) {
      cli::cli_alert_warning(paste(
        "No model passes the adequacy test at alpha = {alpha};",
        "falling back to the global minimum"
      ))
      best_topic <- global_min
      rule <- "global minimum (sequential fallback)"
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
      "Optimal model has {best_topic$topic} topics (selected by {rule})"
    )
  }

  if (do_plot) {
    x_min <- best_topic$topic
    y_min <- best_topic$OpTop
    p1 <- ggplot2::ggplot(Chi_K) +
      ggplot2::geom_line(ggplot2::aes(x = topic, y = OpTop),
                         linewidth = 0.8, color = "royalblue") +
      ggplot2::geom_hline(yintercept = y_min, color = "black", linetype = 2L) +
      ggplot2::geom_vline(xintercept = x_min, color = "black", linetype = 2L) +
      ggplot2::annotate("point", x = x_min, y = y_min,
                        color = "red", shape = 4L, size = 4L) +
      ggplot2::xlab("Topics") + ggplot2::ylab(expression(OpTop[J]^{"K"})) +
      ggplot2::ggtitle("Optimal Topic Plot") +
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
