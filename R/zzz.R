#' @keywords internal
#' @useDynLib OpTop, .registration = TRUE
NULL
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "_OpTop_optimal_topic_core"
  ))
}