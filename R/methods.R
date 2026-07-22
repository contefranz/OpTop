# S3 print and plot methods for the objects OpTop returns. Each print
# method emits a compact, self-explanatory header through cli and returns
# its argument invisibly; the underlying lists are untouched, so
# established dollar-access code keeps working unchanged.

# non-standard evaluation inside ggplot2::aes()
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("K", "metric", "macro", "ci_lo", "ci_hi",
                           "gain", "upper"))
}

#' Printing OpTop objects
#'
#' Compact console summaries for the objects the package returns: models
#' wrapped for OpTop, harmonized partitions, discrepancy-index results,
#' grid tables, held-out results, and moment tests. Each method prints a
#' short header with the quantities a reader checks first and returns the
#' object invisibly; the objects themselves stay plain lists (or a
#' `data.table` for the grid table), so `$`-access and existing code are
#' unaffected.
#'
#' @param x The object to print.
#' @param ... Passed on where a body table is printed (the grid table
#'   method forwards to the `data.table` printer); ignored otherwise.
#'
#' @return `x`, invisibly.
#'
#' @name optop-print
NULL

#' @rdname optop-print
#' @export
print.optop_theta_phi <- function(x, ...) {
  cli::cli_text("{.cls optop_theta_phi}: a topic model for OpTop")
  cli::cli_text("K = {.val {x$K}} topics, {.val {length(x$docs)}} document{?s}, {.val {length(x$terms)}} term{?s}")
  cli::cli_text("theta: {nrow(x$theta)} x {ncol(x$theta)}, phi: {nrow(x$phi)} x {ncol(x$phi)}; rows sum to 1")
  invisible(x)
}

#' @rdname optop-print
#' @export
print.optop_partition <- function(x, ...) {
  J <- length(x$L)
  W <- length(x$vocab)
  total <- x$nonrare_offsets[J + 1]
  rep_ <- x$chisq_min_report
  cli::cli_text("{.cls optop_partition}: harmonized support, {.val {J}} document{?s} x {.val {W}} term{?s} (c = {x$c})")
  cli::cli_text("non-rare entries: {.val {total}} ({sprintf('%.1f', total / max(J, 1))} per document on average)")
  if (rep_$n_excluded > 0) {
    cli::cli_text("Pearson min-bin rule: {.val {rep_$n_excluded}} document{?s} excluded ({sprintf('%.1f', 100 * rep_$share)}% of the corpus)")
  } else {
    cli::cli_text("Pearson min-bin rule: no documents excluded")
  }
  invisible(x)
}

#' @rdname optop-print
#' @export
print.optop_index <- function(x, ...) {
  word <- "r2_word" %in% names(x)
  lvl <- if (word) " at the word level" else ""
  cli::cli_text("{.cls optop_index}: {x$metric} discrepancy index{lvl}, K = {.val {x$K}}")
  if (word) {
    n_valid <- sum(!is.na(x$r2_word))
    cli::cli_text("Micro R^2: {sprintf('%.4f', x$r2_micro_word)} | Macro R^2: {sprintf('%.4f', x$r2_macro_word)}")
    cli::cli_text("V+ = {.val {n_valid}} of {.val {length(x$r2_word)}} word{?s} with a defined index")
  } else {
    macro <- if (!is.null(x$r2_macro)) {
      sprintf(" | Macro R^2: %.4f", x$r2_macro)
    } else {
      ""
    }
    cli::cli_text("Micro R^2: {sprintf('%.4f', x$r2)}{macro}")
    n_valid <- sum(!is.na(x$r2_doc))
    if (x$n_null_excluded > 0) {
      cli::cli_text("J+ = {.val {n_valid}} of {.val {length(x$r2_doc)}} document{?s} ({.val {x$n_null_excluded}} excluded by the null floor)")
    } else {
      cli::cli_text("J+ = {.val {n_valid}} of {.val {length(x$r2_doc)}} document{?s}")
    }
  }
  invisible(x)
}

#' @rdname optop-print
#' @export
print.optop_index_table <- function(x, ...) {
  cli::cli_text("{.cls optop_index_table}: discrepancy indices over the topic grid")
  NextMethod()
  invisible(x)
}

#' @rdname optop-print
#' @export
print.optop_holdout <- function(x, ...) {
  ks <- paste(x$K, collapse = ", ")
  mets <- paste(x$metrics, collapse = ", ")
  J <- nrow(x$scores[[1L]])
  cli::cli_text("{.cls optop_holdout}: held-out discrepancy indices over K = {{{ks}}}")
  cli::cli_text("{.val {J}} evaluation document{?s}, {format(100 * x$conf)}% intervals, metrics: {mets}")
  if (x$stabilize > 0) {
    cli::cli_text("stabilized variant, kappa = {x$stabilize}")
  }
  print(x$summary)
  invisible(x)
}

#' @rdname optop-print
#' @export
print.optop_moment_test <- function(x, ...) {
  lbl <- switch(x$type,
                contrast = "frequency-contrast screen",
                strata = "frequency-strata Wald test",
                fit = "fit-stratified Wald test",
                x$type)
  q <- x$summary$q[1L]
  cli::cli_text("{.cls optop_moment_test}: {lbl} ({.val {q}} instrument{?s})")
  print(x$summary)
  invisible(x)
}

#' Plot a held-out validation result
#'
#' Visualize an [optop_index_holdout()] result. The `"macro"` view draws
#' the held-out Macro index against the number of topics for every
#' evaluated metric, with its confidence band: the working view for
#' locating where the curve flattens. The `"gains"` view draws the
#' adjacent paired gains of one metric with their one-sided upper bounds
#' and the epsilon-adequacy rule of [optop_gain_table()]: the first step
#' out of K whose upper bound falls below `epsilon` selects K, marked by
#' the dashed vertical line.
#'
#' @param x An [optop_index_holdout()] result.
#' @param which `"macro"` (default) or `"gains"`.
#' @param metric For the gains view, the metric to draw (default: the
#'   first evaluated metric); the macro view always draws all of them.
#' @param epsilon,alpha The adequacy tolerance and the one-sided level of
#'   the gains view, passed to [optop_gain_table()].
#' @param ... Ignored.
#'
#' @return The ggplot object, invisibly (the plot is drawn as a side
#'   effect).
#'
#' @examples
#' \donttest{
#' # simulate a corpus with known truth and split train/eval
#' rdirich <- function(n, k) {
#'   g <- matrix(stats::rgamma(n * k, shape = 1), n, k)
#'   g / rowSums(g)
#' }
#' theta <- rdirich(60, 3)
#' phi <- rdirich(3, 120)
#' colnames(phi) <- sprintf("w%03d", 1:120)
#' dtm <- sim_dfm(theta, phi, doc_length = 150, seed = 7)
#'
#' split <- seq_len(40)
#' models <- lapply(2:5, function(k) {
#'   topicmodels::LDA(quanteda::convert(dtm[split, ], to = "topicmodels"),
#'                    k = k, method = "VEM",
#'                    control = list(seed = 500 + k))
#' })
#' baseline <- optop_make_baseline(dtm[split, ])
#' ho <- optop_index_holdout(models, dtm[-split, ], baseline,
#'                           metrics = "deviance")
#' plot(ho)
#' plot(ho, which = "gains", epsilon = 0.01)
#' }
#'
#' @export
plot.optop_holdout <- function(x, which = c("macro", "gains"),
                               metric = NULL, epsilon = 0.01, alpha = 0.05,
                               ...) {
  which <- match.arg(which)
  if (which == "macro") {
    p <- .optop_plot_macro_view(as.data.frame(x$summary), x$conf,
                                "Held-out")
  } else {
    if (is.null(metric)) {
      metric <- x$metrics[1L]
    }
    gt <- optop_gain_table(x, metric = metric, epsilon = epsilon,
                           alpha = alpha)
    p <- .optop_plot_gains_view(as.data.frame(gt$gains), epsilon, alpha,
                                gt$k_hat, metric, "Held-out")
  }
  print(p)
  invisible(p)
}

#' Printing and plotting a cross-fitted result
#'
#' `print.optop_crossfit` summarizes an [optop_crossfit()] run: the
#' configuration, the epsilon-adequacy pick per metric at the default
#' tolerance, and the pooled summary table. `plot.optop_crossfit` draws
#' the same two views as [plot.optop_holdout()], the pooled Macro curve
#' with its fold-replicate confidence band and the cross-fitted
#' adjacent-gains view; `epsilon` and `alpha` recompute the gains at
#' plot time.
#'
#' @param x An [optop_crossfit()] result.
#' @param which `"macro"` (default) or `"gains"`.
#' @param metric For the gains view, the metric to draw (default: the
#'   first evaluated metric).
#' @param epsilon,alpha The adequacy tolerance and one-sided level of
#'   the gains view.
#' @param ... Ignored.
#'
#' @return `x` invisibly for the print method; the ggplot object,
#'   invisibly, for the plot method (drawn as a side effect).
#'
#' @name optop-crossfit-methods
NULL

#' @rdname optop-crossfit-methods
#' @export
print.optop_crossfit <- function(x, ...) {
  ks <- paste(x$K, collapse = ", ")
  mets <- paste(x$metrics, collapse = ", ")
  cli::cli_text("{.cls optop_crossfit}: V-fold cross-fitted discrepancy indices")
  cli::cli_text("{.val {length(x$folds)}} document{?s} in {.val {x$V}} fold{?s}, K = {{{ks}}}, metrics: {mets}")
  cli::cli_text("{format(100 * x$conf)}% intervals from fold-replicate standard errors")
  for (m in x$metrics) {
    kh <- x$gains[[m]]$k_hat
    if (is.finite(kh)) {
      cli::cli_text("epsilon-adequacy pick ({m}, epsilon = {x$gains[[m]]$epsilon}, alpha = {x$gains[[m]]$alpha}): K = {.val {kh}}")
    } else {
      cli::cli_text("epsilon-adequacy pick ({m}): no step qualifies at epsilon = {x$gains[[m]]$epsilon}")
    }
  }
  print(x$summary)
  invisible(x)
}

#' @rdname optop-crossfit-methods
#' @export
plot.optop_crossfit <- function(x, which = c("macro", "gains"),
                                metric = NULL, epsilon = 0.01,
                                alpha = 0.05, ...) {
  which <- match.arg(which)
  if (which == "macro") {
    p <- .optop_plot_macro_view(as.data.frame(x$summary), x$conf,
                                "Cross-fitted")
  } else {
    if (is.null(metric)) {
      metric <- x$metrics[1L]
    }
    if (!metric %in% x$metrics) {
      stop("metric was not evaluated in this cross-fitting run")
    }
    gt <- .optop_crossfit_gains(x, metric, epsilon, alpha)
    p <- .optop_plot_gains_view(as.data.frame(gt$gains), epsilon, alpha,
                                gt$k_hat, metric, "Cross-fitted")
  }
  print(p)
  invisible(p)
}

# Shared view builders of the holdout and crossfit plot methods: one
# Macro-curve view over the summary table, one adjacent-gains view over a
# gains table carrying K, gain, and upper.
.optop_plot_macro_view <- function(s, conf, prefix) {
  ggplot2::ggplot(
    s, ggplot2::aes(x = K, y = macro, color = metric, fill = metric)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      alpha = 0.15, color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::labs(
      x = "number of topics K",
      y = sprintf("%s Macro index", tolower(prefix)),
      color = "metric", fill = "metric",
      title = sprintf("%s Macro discrepancy index", prefix),
      subtitle = sprintf("%g%% pointwise intervals", 100 * conf)
    )
}

.optop_plot_gains_view <- function(g, epsilon, alpha, k_hat, metric,
                                   prefix) {
  p <- ggplot2::ggplot(g, ggplot2::aes(x = K, y = gain)) +
    ggplot2::geom_hline(yintercept = epsilon, linetype = 3) +
    ggplot2::geom_hline(yintercept = 0, linetype = 1,
                        color = "grey70", linewidth = 0.3) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = K, y = gain, yend = upper),
      linewidth = 0.4
    ) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::labs(
      x = "step out of K (to the next point of the grid)",
      y = sprintf("%s adjacent gain", tolower(prefix)),
      title = sprintf("Adjacent gains, %s index", metric),
      subtitle = sprintf(
        "one-sided %g%% upper bounds; dotted line: epsilon = %g",
        100 * (1 - alpha), epsilon
      )
    )
  if (is.finite(k_hat)) {
    p <- p + ggplot2::geom_vline(xintercept = k_hat, linetype = 2)
  }
  p
}
