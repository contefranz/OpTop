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
#'                               optimal_model = test1,
#'                               best_match = BestMatch,
#'                               least_match = LeastMatch )
#' }
#' @seealso \code{\link[topicmodels]{LDA}} \code{\link[data.table]{data.table}}
#' @references Lewis, C. and Grossetti, F. (2019 - forthcoming):\cr
#' A Statistical Approach for Optimal Topic Model Identification.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}
#' @author Craig M. Lewis \email{craig.lewis@@owen.vanderbilt.edu}
#' @import data.table
#' @export

agg_topic_stability <- function( lda_models, optimal_model, 
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
  if ( ( !is.data.frame( optimal_model ) || !is.data.table( optimal_model ) ) &&
       !is.numeric( optimal_model ) && !is.LDA_VEM( optimal_model ) ) {
    stop( paste( "optimal_model must be either 1. an integer identifying",
                 "the number of topics which best fits the corpus",
                 "2. a data.table/data.frame as returned by .optimal_model()",
                 "3. an LDA_VEM obect as computed by topicmodels::LDA()" ) )
  }
  if ( !is.numeric( threshold ) ) {
    stop( "threshold must be either a numeric or an integer" )
  }
  if ( !is.logical( convert ) ) {
    stop( "convert must be either TRUE or FALSE" )
  }
  if ( is.numeric( optimal_model ) ) {
    .optimal_model <- optimal_model
  } else if ( is.data.table( optimal_model ) || is.data.frame( optimal_model ) ) {
    cat( "optimal_model is a data.table or a data.frame.",
         "Extracting information about optimal model...\n" )
    .optimal_model <- optimal_model[ which.min( chi_std ), topic ]
  } else if ( is.LDA_VEM( optimal_model ) ) {
    cat( "optimal_model is a LDA_VEM object.", 
         "Extracting information about the optimal model...\n" )
    dtw_best <- optimal_model@gamma
    tww_best <- t( exp( optimal_model@beta ) )
    .optimal_model <- ncol( dtw_best )
  }
  cat( "best model has", .optimal_model, "topics\n" )
  
  tic <- proc.time()
  
  k_end <- max( sapply( lda_models, function( x ) x@k ) )
  best_pos <- which( sapply( lda_models, function( x ) x@k ) == .optimal_model )
  
  if ( length( best_pos ) == 0 ) {
    stop( paste( "There is no optimal model in lda_models.",
                 "This could be either due to a wrong specification of",
                 "argument optimal_model or",
                 "if optimal_model is a data.table, the optimal model cannot be found",
                 "in the list lda_models." ) )
  }
  if ( !is.LDA_VEM( optimal_model ) ) {
    # extracting information from best model
    dtw_best <- lda_models[[ best_pos ]]@gamma
    tww_best <- t( exp( lda_models[[ best_pos ]]@beta ) )
  }
  n_doc <- nrow( dtw_best )
  loop_sequence <- (best_pos + 1L):length( lda_models )  
  k = 1L
  out <- data.table()
  cat( "---\n" )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  for ( i_mod in loop_sequence ) {
    i_pos <- i_mod - .optimal_model + 1L
    dtw <- lda_models[[ i_mod ]]@gamma
    tww <- t( exp( lda_models[[ i_mod ]]@beta ) )
    current_k <- ncol( tww )
    n_tww <- nrow( tww )
    p_tww <- ncol( tww )
    cat( "---\n" )
    cat( "# # # Processing LDA with k =", current_k, "\n" )
    
    i <- current_k - .optimal_model
    MostSim <- intersect( 1L:current_k, best_match[ i, ] )
    MostSimDTW <- dtw[ , MostSim ]
    MostSimTWW <- tww[ , MostSim ]
    
    out_in <- data.table()
    cat( "--> Processing documents\n" )
    for ( j_doc in 1L:n_doc ) {
      DTW_best <- matrix( dtw_best[ j_doc, ], 
                          nrow = n_tww,
                          ncol = .optimal_model )
      DTW_most <- matrix( MostSimDTW[ k, ], 
                          nrow = n_tww,
                          ncol = ncol( MostSimDTW ) )
      
      tww_dtw_best <- DTW_best * tww_best
      ones_tww_dtw_best <- matrix( 1., nrow = ncol( tww_dtw_best ), ncol = 1L )
      # totprob_best <- sum( tww_dtw_best %*% rep( 1L, ncol( tww_dtw_best ) ) )
      totprob_best <- sum( tww_dtw_best %*% ones_tww_dtw_best, na.rm = TRUE )
      tww_dtw_best2 <- tww_dtw_best %*% rep( 1L, ncol( tww_dtw_best ) ) / totprob_best
      
      tww_dtw <- DTW_most * MostSimTWW
      totprob <- sum( tww_dtw %*% rep( 1L, ncol( tww_dtw ) ) )
      tww_dtw2 <- tww_dtw %*% rep( 1L, ncol( tww_dtw ) ) / totprob
      
      X <- cbind( tww_dtw_best2, tww_dtw2 )
      BestPair <- apply( X, 2L, function( x ) sort( x, decreasing = TRUE ) )
      icut <- base::which.min( abs( BestPair[ , 1L ] - threshold ) )
      if ( icut > 250 ) {
        sum_overbest <- apply( BestPair[ (icut + 1L):nrow( BestPair ), ], 2L, sum )
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