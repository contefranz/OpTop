#' OpTop: detect the optimal number of topics from a pool of LDA models
#' 
#' OpTop is an \code{R} package that implements the testing approach described in 
#' the paper A Statistical Approach for Optimal Topic Model Identification
#' by Lewis and Grossetti (2019). 
#' 
#' OpTop introduces a set of parametric tests to identify the optimal number 
#' of topics in a collection of LDA models. OpTop also includes several 
#' tests to explore topic stability and redundancy.
#'
#' @docType package
#' @name OpTop
#' @useDynLib OpTop, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
#> NULL
