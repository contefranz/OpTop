# if ( getRversion() >= "2.15.1" ) {
#   utils::globalVariables( c( "id_word", "id_doc", "word_proportions", 
#                              "word_prop",
#                              "word_sum", "check", "topics",
#                              ".", "chisquare", "chisquare_mod",
#                              "row_cut", "chi_sum", "word_prop_hat",
#                              "word_prop_hat_cum", "pval" ) )
# }
#' Simulate a document-feature-matrix from a LDA specification
#'
#' Simulate a document-feature-matrix from estimated DTW and TWW. The number of topics as well as 
#' the number of documents are inferred from the LDA parameters.
#'
#' @param DTW A matrix or data.frame with Document-Topic-Weights.
#' @param TWW A matrix or data.frame with Topic-Word-Weights.
#' @param doc_lenght A vector containing the desired document length as total number of word counts.
#' @param alpha Parameter of the Dirichlet distribution for topics over documents.
#' @param seed Input to \code{set.seed}.
#' 
#' @return A \code{\link[quanteda]{dfm}} object.
#' 
#' @seealso \code{\link[topicmodels]{LDA}} \code{\link[quanteda]{dfm}}
#' @references Lewis, C. and Grossetti, F. (2019 - forthcoming):\cr
#' A Statistical Approach for Optimal Topic Model Identification.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}
#' @author Craig M. Lewis \email{craig.lewis@@owen.vanderbilt.edu}
#' @import data.table ggplot2
#' @importFrom LDATS sim_LDA_data
#' @importFrom quanteda as.dfm
#' @export

sim_dfm <- function( DTW, TWW, doc_lenght, alpha = NULL, seed = NULL) {
  
  if ( all( !is.matrix( DTW ), !is.data.frame( DTW ) ) ) {
    stop( "DTW must be either a matrix or a data.frame" )
  }
  if ( all( !is.matrix( TWW ), !is.data.frame( TWW ) ) ) {
    stop( "TWW must be either a matrix or a data.frame" )
  }
  if ( all( !is.integer( doc_lenght ), !is.numeric( doc_lenght ) ) ) {
    stop( "doc_lenght must be either an integer or a numeric vector" )
  }
  
  out <- sim_LDA_data( N = doc_lenght,
                       Beta = TWW,
                       alpha = alpha,
                       Theta = DTW,
                       seed = seed )
  
  outdfm <- as.dfm( out )
  return( outdfm )
  
}

