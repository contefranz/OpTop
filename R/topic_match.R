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
  if ( ( !is.data.frame( optimal_model ) || !is.data.table( optimal_model ) ) &&
       !is.numeric( optimal_model ) && !is.LDA_VEM( optimal_model ) ) {
    stop( paste( "optimal_model must be either 1. an integer identifying",
                 "the number of topics which best fits the corpus",
                 "2. a data.table/data.frame as returned by optimal_topic()",
                 "3. an LDA_VEM obect as computed by topicmodels::LDA()" ) )
  }
  if ( is.numeric( optimal_model ) ) {
    .optimal_model <- optimal_model
  } else if ( is.data.table( optimal_model ) || is.data.frame( optimal_model ) ) {
    cat( "optimal_model is a data.table or a data.frame.",
         "Extracting information about optimal model...\n" )
    .optimal_model <- optimal_model[ which.min( OpTop ), topic ]
  } else if ( is.LDA_VEM( optimal_model ) ) {
    cat( "optimal_model is a LDA_VEM object.", 
         "Extracting information about the optimal model...\n" )
    dtw_best <- optimal_model@gamma
    tww_best <- t( exp( optimal_model@beta ) )
    .optimal_model <- ncol( dtw_best )
  }
  
  if ( .optimal_model == lda_models[[ length(lda_models ) ]]@k ) {
    message("Optimal model is already the last one in lda_models. There is nothing to compute above that.")
    return( NULL )
  }
  cat( "best model has", .optimal_model, "topics\n" )  
  tic <- proc.time()
  
  k_end <- max( sapply( lda_models, function( x ) x@k ) )
  # find the element corresponding to the best topic
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
    tww_best <- t( exp( lda_models[[ best_pos ]]@beta ) )
  }
  # Normalizing by scaling by vector norms
  tww_best_norm <- .Call(`_OpTop_normalize_columns`, tww_best)
  
  loop_sequence <- (best_pos + 1L):length( lda_models )
  min_loop <- min( loop_sequence )
  l_loop <- length( loop_sequence )
  BestMatch  <- LeastMatch <- vector( "list", length = l_loop )
  # initialize topic position for BestMatch and LeastMatch
  k = 1L
  cat( "---\n" )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  for ( i_mod in loop_sequence ) {
    # i_pos <- i_mod - .optimal_model + 1L
    i_pos <- i_mod - min_loop + 1L
    
    # reading tww
    tww <- t( exp( lda_models[[ i_mod ]]@beta ) )
    current_k <- ncol( tww )
    # n_tww <- nrow( tww )
    # p_tww <- ncol( tww )
    cat( "---\n" )
    cat( "# # # Processing LDA with k =", current_k, "\n" )
    # Normalizing by scaling by vector norms
    tww_norm <- .Call(`_OpTop_normalize_columns`, tww)
    
    # this is the matrix multiplication to get the Cosine Similarities
    # between the best set of topics and the targe ones
    CosSim <- t( tww_best_norm ) %*% tww_norm
    maxval <- apply( CosSim, 1L, max )
    maxsim <- apply( CosSim, 1L, which.max )
    
    # vectorize CosSim matrix
    CosSimVec <- as.numeric( CosSim )
    # Identify minimum cosine similarity threshhold.
    # If CosSim exceeds "thresh" it is too similar 
    # to be treated as a distinct factor
    # Matlab uses a weighting of 0.70 which is difficult to replicate here
    if ( var_correction ) {
      thresh <- min( mean( CosSimVec ) + 
                       2.58*stats::sd( CosSimVec, na.rm = TRUE ) )
    } else {
      n = length( CosSimVec )
      thresh <- min( mean( CosSimVec ) + 
                       2.58*stats::sd( CosSimVec, na.rm = TRUE ) * 
                       sqrt((n - 1L)/n), na.rm = TRUE)
    }
    
    tops <- 1L:current_k
    # pre-allocate object
    above_thresh <- matrix( as.integer( CosSim > thresh ), 
                            nrow = .optimal_model, 
                            ncol = current_k )
    top_mat      <- matrix( rep( 1L:current_k, each = .optimal_model ), 
                            nrow = .optimal_model, 
                            ncol = current_k )
    
    # Identify factors that are highly similar to at least one base model factor
    CosSimCheck <- above_thresh * top_mat
    # A single factor may be highly similar to multiple base model 
    # factors. This step creates a vector of unique factor identifiers
    BigCosSim <- intersect( tops, CosSimCheck )
    # Merge all factors considered to be most similar
    Inform   <- sort( union( maxsim, BigCosSim ) )
    # Identify all factors that are considered to be "uninformative"
    Uninform <- sort( setdiff( tops, Inform ) )
    
    x <- length( Inform ) + 1L
    BestMatch[[ i_pos ]] <- c( current_k, Inform )
    x <- length( Uninform ) + 1L
    if ( x > 1L ) {
      LeastMatch[[ i_pos ]] <- c( current_k, Uninform )
    }
  } 
  
  # converting lists to matrix accounting accounting for different
  # vectors lengths
  best_max_length  <- max( sapply( BestMatch, length ), na.rm = TRUE ) 
  least_max_length <- max( sapply( LeastMatch, length ), na.rm = TRUE ) 
  
  for ( i_maxing in seq_along( BestMatch ) ) {
    length( BestMatch[[ i_maxing ]] ) <- best_max_length
    length( LeastMatch[[ i_maxing ]] ) <- least_max_length
  }
  
  BestMatch_out <- matrix( unlist( BestMatch ), 
                           nrow = length( BestMatch ), 
                           ncol = best_max_length, byrow = TRUE )
  LeastMatch_out <- matrix( unlist( LeastMatch ), 
                            nrow = length( LeastMatch ), 
                            ncol = least_max_length, byrow = TRUE )
  
  cat( "---\n" )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  
  toc <- proc.time()
  runtime <- toc - tic
  cat( "---\n" )
  cat( "Function took:", runtime[ 3L ], "sec.\n" )
  cat( "---\n" )
  return( list( BestMatch = BestMatch_out, LeastMatch = LeastMatch_out ) )
  
}
