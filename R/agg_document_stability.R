if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "chi_sq_std", "df", ".") )
}
#' Compute aggregate document stability
#'
#' Detects informative and uninformative components to compute aggregate document
#' stability.
#' 
#' @inheritParams optimal_topic
#' @inheritParams topic_stability
#' @inheritParams agg_topic_stability
#' @examples
#'\dontrun{
#' test4 <- agg_topic_stability( lda_models = lda_list,
#'                               best_k = test1,
#'                               best_match = BestMatch,
#'                               least_match = LeastMatch )
#' }
#' @seealso \code{\link[topicmodels]{LDA}} \code{\link[data.table]{data.table}}
#' @references Lewis, C. and Grossetti, F. (2019 - forthcoming):\cr
#' A Statistical Approach for Optimal Topic Model Identification.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}.
#' @import data.table
#' @export

agg_document_stability <- function( lda_models, word_proportions, 
                                    optimal_model, 
                                    threshold = 0.00075, convert = FALSE ) {
  
  if ( !is.list( lda_models ) ) {
    stop( "lda_models must be a list" )
  }
  if ( length( lda_models ) == 1L ){
    stop( paste( "length(lda_models) = 1.",
                 "This is strange since the test should be perfomed",
                 "on multiple LDA models." ) )
  }
  if ( !all( sapply( lda_models, is.LDA_VEM ) ) ) {
    stop( paste( "lda_models must contain LDA_VEM obects as computed",
                 "by topicmodels::LDA()" ) )
  }
  if( !is.data.table( word_proportions ) ) {
    stop( "word_proportions must be a data.table" )
  }
  if ( !is.numeric( threshold ) ) {
    stop( "threshold must be either a numeric or an integer" )
  }
  if ( !is.logical( convert ) ) {
    stop( "convert must be either TRUE or FALSE" )
  }
  
  tic <- proc.time()
  # compute the size of vocabulary detected in each document as:
  size_corpus <- nrow( word_proportions )
  size_vocabulary <- nrow( word_proportions[ , .N, by = id_word ] )
  n_docs <- nrow( word_proportions[ , .N, by = id_doc ] )
  
  
  
  
}