if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "chi_std", "topic", "pchisq" ) )
}
#' Compute topic stability for over-optimal topic specifications
#'
#' Implements fast chi-square like tests to evaluate the stability of redundant
#' topics. Both point-wise chi-square statistics and aggregated one can be 
#' selected.
#' 
#' @inheritParams optimal_topic
#' @param best_k The optical topic model. It can be either an integer giving
#' the number of topics for the optimal topic model or a
#' \code{\link[data.table]{data.table}}/\code{\link[base]{data.frame}} as
#' returned by \code{\link[OpTop]{optimal_topic}}. In the latter case, the
#' function automatically finds the optimal model.
#' @param alpha Alpha level to identify informative words from the Cumulative
#' Distribution Function over the cosine similarities in the Topic Word Weights
#' matrix. Default to 0.05.
#' @param test A character specifying what the function returns. By default,
#' \code{"aggregated"} is chosen so that the aggregated chi-square statistics 
#' as specified in Test 3 of the paper is returned. Other specifications include
#' \code{"single"} for point-wise chi-square statistics as given by Test 2 and 
#' \code{"both"} to get both Tests 2 and 3.
#' See 'Details' for more information.
#' @param compute_res Determines whether or not the function has to compute
#' cosine similarities for unmatched topic.
#' @details This function implements both Tests 2 and 3 as defined in 
#' Lewis and Grossetti (2019). Test 2 evaluates the point-wise chi-square statistics
#' for every over-optimal topic specifications. Test 3 evaluates the aggregated 
#' stability of over-optimal topic specifications by summing each 
#' point-wise contribution. See 'Value' to understand how \code{topic_stability} 
#' returns the results.
#' @return Either a data.table, a matrix or a list as specified by argument 
#' 
#' If \code{test = "aggregated"}, a data.table representing Test 3 
#' is returned with the following columns:
#' \item{\code{dfs}}{An integer giving the degrees of freedom.}
#' \item{\code{chisq_std}}{A numeric giving the standardized chi-square statistics.}
#' 
#' ---
#' 
#' If \code{test = "single"}, a matrix of size 
#' \code{[max{k} - best_k] X [best_k]} representing Test 2 is returned.
#' 
#' ---
#' 
#' If \code{test = "both"}, a named list with both Tests 2 and 3 is returned.
#' @examples
#' \dontrun{
#' test2 <- topic_stability( lda_models = lda_list,
#'                           best_k = test1,
#'                           threshold = 0.00075,
#'                           alpha = 0.05 )
#' }
#' @seealso \code{\link[topicmodels]{LDA}} \code{\link[data.table]{data.table}}
#' @references Lewis, C. and Grossetti, F. (2019 - forthcoming):\cr
#' A Statistical Approach for Optimal Topic Model Identification.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}.
#' @import data.table
#' @importFrom stats pchisq
#' @export

topic_stability <- function( lda_models, best_k,
                             threshold = 0.00075, alpha = 0.05,
                             test = c( "aggregated", "single", "both" ),
                             compute_res = FALSE, convert = FALSE ) {

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
  if ( !is.numeric( threshold ) ) {
    stop( "threshold must be either a numeric or an integer" )
  }
  if ( !is.numeric( alpha ) ) {
    stop( "alpha must be a numeric" )
  }
  if ( !is.character( test ) ) {
    stop( paste("test must be a character string specifying either 'aggregated'",
                "'single', or 'both' to choose what the topic_stability() returns"))
  }
  if ( compute_res ) {
    warning( "Argument compute_res is deprecated and will be removed soon!" )
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
  tww_best <- t( exp( lda_models[[ best_pos ]]@beta ) )
  n_best <- nrow( tww_best )
  p_best <- ncol( tww_best )
  if ( p_best != best_k ) {
    stop( "Wrong identification of optimal topic model!" )
  }

  # Normalizing by scaling by vector norms
  tww_best_norm <- norm_tww( tww_best )
  
  # pre-allocating required objects
  BestTopics    <- matrix( 0, nrow = (k_end - best_k), ncol = best_k )
  Chi_k         <- matrix( 0, nrow = (k_end - best_k), ncol = (p_best + 1L) )
  Chi_k[ , 1L ] <- (best_k + 1L):k_end
  Chi_K         <- matrix( 0, nrow = (k_end - best_k), ncol = 3L )
  Chi_K[ , 1L ] <- 1L:nrow( Chi_K )
  p             <- matrix( 0, nrow = (k_end - best_k), ncol = (p_best + 1L) )
  CosSimHold    <- NULL

  # this loops over the list of models which have to be sorted
  # the loops starts right after the best model
  loop_sequence <- (best_pos + 1L):length( lda_models )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  for ( i_mod in loop_sequence ) {
    i_pos <- i_mod - best_k + 1L
    # reading tww
    tww <- t( exp( lda_models[[ i_mod ]]@beta ) )
    current_k <- ncol( tww )
    n_tww <- nrow( tww )
    p_tww <- ncol( tww )
    cat( "---\n" )
    cat( "# # # Processing LDA with k =", current_k, "\n" )

    # Normalizing by scaling by vector norms
    tww_norm <- norm_tww( tww )
    
    # this is the matrix multiplication to get the Cosine Similarities
    # between the best set of topics and the targe ones
    CosSim <- t( tww_best_norm ) %*% tww_norm
    # transpose the result so that N reflects the number of topics under
    # investigation
    minCosSim <- t( CosSim )
    CumDF <- 0

    # and here begins the second loop
    cat( "--> Finding topics with highest cosine similarity\n" )
    for ( j in 1L:p_best ) {
      # find the max cosine similarity and return the value and its index
      maxval <- max( CosSim[ j, ] )
      maxsim <- which.max( CosSim[ j, ] )
      BestTopics[ (current_k - best_k), j ] <- maxsim
      BestPair <- matrix( data = c( tww_best[ , j ],
                                    tww[ , maxsim ] ),
                          ncol = 2L )
      # sot by descending order each column of BestPair
      BestPair <- apply( BestPair, 2L, function( x ) sort( x, decreasing = TRUE ) )
      CumProb <- 1 - cumsum( BestPair[ , 1L ] )
      # R cannot manage extemely small numbers. We have to impute 1
      # as first element of CumProb and kill last element
      CumProb <- c( 1, CumProb[ -length( CumProb ) ] )
      BestPair <- cbind( BestPair, unname( CumProb ) )

      icut <- base::which.min( abs( BestPair[ , 1L ] - threshold ) )
      ibin <- base::which.min( abs( BestPair[ , 3L ] - alpha ) )
      ibestcut <- min( icut, ibin )
      sum_overbest <- apply( BestPair[ (ibestcut + 1L):n_tww, ], 2L, sum )
      AggBestPair <- rbind( BestPair[ 1L:ibestcut, ], unname( sum_overbest ) )
      numerator <- ( AggBestPair[ , 1L ] - AggBestPair[ , 2L ] )^2L
      denominator <- AggBestPair[ , 1L ]
      Chi_k[ i_pos, j + 1L ] <- ibestcut * sum( numerator / denominator )
      p[ i_pos, j ] <- 1 - pchisq( Chi_k[ i_pos, j + 1L ], ibestcut )

    }
    if ( compute_res ) {
      ResTopics <- setdiff( 1L:current_k, BestTopics[ (current_k - best_k), ] )
      for ( i_res in seq_along( ResTopics ) ) {
        TopicsMatch <- CosSim[ , ResTopics[ i_res ] ]
        CosSimHold <- cbind( CosSimHold, current_k - 1L,
                             ResTopics[ i_res ], t( TopicsMatch ) )
      }
    }
  }

  Chi_K[ , 2L ] <- apply( Chi_k, 1L, sum ) / best_k
  # not sure if we want to divide by 1:108 or by the real k given
  # by topics over optimal
  Chi_K[ , 3L ] <- Chi_K[ , 2L ] / Chi_K[ , 1L ]
  colnames( Chi_K ) <- c( "dfs", "chisq", "chisq_std" )
  
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  if ( convert ) {
    Chi_K <- as.data.table( Chi_K )
  } else {
    Chi_K <- as.data.frame( Chi_K )
  }
  
  test = match.arg( test )
  if ( test == "aggregated" ) {
    toc <- proc.time()
    runtime <- toc - tic
    cat( "---\n" )
    cat( "Function took:", runtime[ 3L ], "sec.\n" )
    cat( "---\n" )
    return( Chi_K )
  } else if ( test == "single" ) {
    toc <- proc.time()
    runtime <- toc - tic
    cat( "---\n" )
    cat( "Function took:", runtime[ 3L ], "sec.\n" )
    cat( "---\n" )
    return( Chi_k )
  } else if ( test == "both" ) {
    toc <- proc.time()
    runtime <- toc - tic
    cat( "---\n" )
    cat( "Function took:", runtime[ 3L ], "sec.\n" )
    cat( "---\n" )
    return_list <- list( test2 = Chi_k, test3 = Chi_K )
    return( return_list )
  }
}
