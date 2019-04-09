if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "chi_sq_std", "df", ".") )
}
#' Compute aggregate topic stability
#'
#' Detects informative and uninformative components to compute aggregate topic
#' stability.
#' 
#' @inheritParams optimal_topic
#' @inheritParams topic_stability
#' @param best_match A data.table as computed by \code{\link[OpTop]{topic_match}}.
#' This contains the models with highest cosine similarity.
#' @param least_match A data.table as computed by \code{\link[OpTop]{topic_match}}.
#' This contains the models with lowest cosine similarity.
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

agg_topic_stability <- function( lda_models, best_k, 
                                 best_match, least_match,
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
  if ( ( !is.data.frame( best_k ) || !is.data.table( best_k ) ) &&
       !is.numeric( best_k ) ) {
    stop( paste( "best_k must be either an integer identifying",
                 "the number of topics which best fits the corpus",
                 "or a data.table/data.frame as returned by",
                 "optimal_topic()" ) )
  }
  if ( !is.numeric( best_k ) ) {
    cat( "best_k is a data.table or a data.frame. Extracting best model...\n" )
    best_k <- best_k[ which.min( chi_std ), topic ]
    cat( "best model has", best_k, "topics\n" )
  }
  
  tic <- proc.time()
  
  n_topics <- length( lda_models )
  k_end <- max( sapply( lda_models, function( x ) x@k ) )
  
  # find the element corresponding to the best topic
  best_pos <- which( sapply( lda_models, function( x ) x@k ) == best_k )
  if ( length( best_pos ) == 0 ) {
    stop( paste( "There is no optimal model in lda_models.",
                 "This could be either due to a wrong specification of",
                 "argument best_k or",
                 "if best_k is a data.table, the optimal model cannot be found",
                 "in the list lda_models." ) )
  }
  # extracting information from best model
  dww_best <- lda_models[[ best_pos ]]@gamma
  n_doc <- nrow( dww_best )
  tww_best <- t( exp( lda_models[[ best_pos ]]@beta ) )
  n_best <- nrow( tww_best )
  p_best <- ncol( tww_best )
  if ( p_best != best_k ) {
    stop( "Wrong identification of optimal topic model!" )
  }
  
  loop_sequence <- (best_pos + 1L):length( lda_models )
  l_loop <- length( loop_sequence )
  k = 1L
  out <- data.table()
  cat( "---\n" )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  for ( i_mod in loop_sequence ) {
    i_pos <- i_mod - best_k + 1L
    
    dww <- lda_models[[ i_mod ]]@gamma
    tww <- t( exp( lda_models[[ i_mod ]]@beta ) )
    current_k <- ncol( tww )
    n_tww <- nrow( tww )
    p_tww <- ncol( tww )
    cat( "---\n" )
    cat( "# # # Processing LDA with k =", current_k, "\n" )
    
    i <- current_k - best_k
    MostSim <- intersect( 1L:current_k, best_match[ i, ] )
    MostSimDWW <- dww[ , MostSim ]
    MostSimTWW <- tww[ , MostSim ]
    
    out_in <- data.table()
    cat( "--> Processing documents\n" )
    for ( j_doc in 1L:n_doc ) {
      DWW_best <- matrix( dww_best[ j_doc, ], 
                          nrow = n_tww,
                          ncol = p_best )
      DWW_most <- matrix( MostSimDWW[ k, ], 
                          nrow = n_tww,
                          ncol = ncol( MostSimDWW ) )
      tww_dww_best <- DWW_best * tww_best
      totprob_best <- sum( tww_dww_best %*% rep( 1L, ncol( tww_dww_best ) ) )
      tww_dww_best2 <- tww_dww_best %*% rep( 1L, ncol( tww_dww_best ) ) / totprob_best
      
      tww_dww <- DWW_most * MostSimTWW
      totprob <- sum( tww_dww %*% rep( 1L, ncol( tww_dww ) ) )
      tww_dww2 <- tww_dww %*% rep( 1L, ncol( tww_dww ) ) / totprob
      
      X <- cbind( tww_dww_best2, tww_dww2 )
      BestPair <- apply( X, 2L, function( x ) sort( x, decreasing = TRUE ) )
      icut <- base::which.min( abs( BestPair[ , 1L ] - threshold ) )
      if ( icut > 250 ) {
        sum_overbest <- apply( BestPair[ (icut + 1L):n_doc, ], 2L, sum )
        AggBestPair <- rbind( BestPair[ 1L:icut, ], unname( sum_overbest ) )
        
        numerator <- ( AggBestPair[ , 1L ] - AggBestPair[ , 2L ] )^2L
        denominator <- AggBestPair[ , 1L ]
        chi_sq <- icut * sum( numerator / denominator, na.rm = TRUE )
        
        out_doc <- data.table( topic = current_k,
                               id_doc = j_doc,
                               df = icut,
                               chi_sq = chi_sq )
        out_in <- rbindlist( list( out_in, out_doc ) )
      }
    }
    out <- rbindlist( list( out, out_in ) )
  }
  out[ , chi_sq_std := chi_sq / df ]
  
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  if ( convert ) {
    cat( "Converting to data.frame\n" )
    setDF( out )
  }
  toc <- proc.time()
  runtime <- toc - tic
  cat( "---\n" )
  cat( "Function took:", runtime[ 3L ], "sec.\n" )
  cat( "---\n" )
  return( out[] )
}