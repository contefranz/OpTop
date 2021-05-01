if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "chi_sq_std", "df" ) )
}
#' Define informative and uninformative components
#'
#' For above-optimal topics, detect and extract informative and uninformative
#' components in terms of cosine similarities with the optimal topic model.
#' 
#' @inheritParams optimal_topic
#' @inheritParams topic_stability
#' @param var_correction Use the unbiased estimator of the co-variance for 
#' i.i.d. observations by applying n - 1 in the denominator. Default is \code{TRUE}.
#' 
#' @return A named list with the informative and uninformative components given
#' as matrices.
#' 
#' @examples
#' \dontrun{
#' out3 <- topic_match( lda_models = lda_list, 
#'                      optimal_model = test1, 
#'                      var_correction = TRUE )
#' }
#' @seealso \code{\link[topicmodels]{LDA}} \code{\link[data.table]{data.table}}
#' @references Lewis, C. and Grossetti, F. (2019 - forthcoming):\cr
#' A Statistical Approach for Optimal Topic Model Identification.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}
#' @author Craig M. Lewis \email{craig.lewis@@owen.vanderbilt.edu}
#' @import data.table
#' @export


topic_match <- function( lda_models, optimal_model, var_correction = TRUE ) {
  
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
  if ( !is.numeric( optimal_model ) || optimal_model < 2 ) {
    stop("optimal_model must be a number greater than 2")
  }
  if ( optimal_model == lda_models[[ length(lda_models ) ]]@k ) {
    message("Optimal model is already the last one in lda_models. There is nothing to compute above that.")
    return( NULL )
  }
  
  tic <- proc.time()
  
  # find the element corresponding to the best topic
  best_pos <- which( sapply( lda_models, function( x ) x@k ) == optimal_model )
  # check that optimal_models does correspond to a real LDA model in lda_models
  if ( length( best_pos ) == 0 ) {
    stop("optimal_model does not correspond to any topic number in lda_models")
  }
  
  cat( "---\n" )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  out = .Call(`_OpTop_topic_match_core`, lda_models, best_pos, optimal_model, var_correction)
  cat( "---\n" )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  
  toc <- proc.time()
  runtime <- toc - tic
  cat( "---\n" )
  cat( "Function took:", runtime[ 3L ], "sec.\n" )
  cat( "---\n" )
  return(out)
  
}
