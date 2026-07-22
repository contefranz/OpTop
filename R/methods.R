# S3 print methods for the objects OpTop returns. Each method prints a
# compact, self-explanatory header through cli and returns its argument
# invisibly; the underlying lists are untouched, so established
# dollar-access code keeps working unchanged.

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
