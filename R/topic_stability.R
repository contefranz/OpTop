if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "chi_std", "topic", "pchisq", "dfs", "chisq_std" ) )
}
#' Compute topic stability for over-optimal topic specifications
#'
#' Implements fast chi-square like test to evaluate the stability of redundant
#' topics.
#' 
#' @inheritParams optimal_topic
#' @param optimal_model The optimal topic model. It can be either an integer giving
#' the number of topics for the optimal topic model or a
#' \code{\link[data.table]{data.table}}/\code{\link[base]{data.frame}} as
#' returned by \code{\link[OpTop]{optimal_topic}}. In the latter case, the
#' function automatically finds the optimal model. This argument can also be
#' a \code{\link[topicmodels]{LDA-class}} object of type \code{LDA_VEM} as 
#' returned by \code{\link[topicmodels]{LDA}}.
#' @param alpha Alpha level to identify informative words from the Cumulative
#' Distribution Function over the cosine similarities in the Topic Word Weights
#' matrix. Default to 0.05.
#' @details This function implements Test 3 as defined in 
#' Lewis and Grossetti (2019). Test 3 evaluates the aggregated 
#' stability of over-optimal topic specifications by summing each 
#' point-wise contribution. See 'Value' to understand how \code{topic_stability} 
#' returns the results.
#' @return A \code{data.table} containing the following columns:
#'
#' \item{\code{topic}}{An integer giving the number of topics.}
#' \item{\code{dfs}}{An integer giving the degrees of freedom.}
#' \item{\code{chisq}}{A numeric giving the chi-square statistic.}
#' \item{\code{chisq_std}}{A numeric giving the standardized chi-square statistic.}
#' @examples
#' \dontrun{
#' test2 <- topic_stability( lda_models = lda_list,
#'                           optimal_model = out,
#'                           threshold = 0.00075,
#'                           alpha = 0.05 )
#' }
#' @seealso \code{\link[topicmodels]{LDA}} \code{\link[data.table]{data.table}}
#' @references Lewis, C. and Grossetti, F. (2019 - forthcoming):\cr
#' A Statistical Approach for Optimal Topic Model Identification.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}
#' @author Craig M. Lewis \email{craig.lewis@@owen.vanderbilt.edu}
#' @import data.table ggplot2
#' @importFrom stats pchisq
#' @export

topic_stability <- function( lda_models, optimal_model,
                             threshold = 0.00075, alpha = 0.05,
                             do_plot = TRUE, convert = FALSE ) {
  
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
  if ( !is.numeric( threshold ) ) {
    stop( "threshold must be either a numeric or an integer" )
  }
  if ( !is.numeric( alpha ) ) {
    stop( "alpha must be a numeric" )
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
  
  # find the element corresponding to the best topic
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
  tww_best_norm <- norm_tww( tww_best )
  
  # pre-allocating required objects
  BestTopics    <- matrix( 0, nrow = (k_end - .optimal_model), ncol = .optimal_model )
  Chi_k         <- matrix( 0, nrow = (k_end - .optimal_model), ncol = (.optimal_model + 1L) )
  Chi_k[ , 1L ] <- (.optimal_model + 1L):k_end
  Chi_K         <- matrix( 0, nrow = (k_end - .optimal_model), ncol = 3L )
  Chi_K[ , 1L ] <- 1L:nrow( Chi_K )
  p             <- matrix( 0, nrow = (k_end - .optimal_model), ncol = (.optimal_model + 1L) )
  CosSimHold    <- NULL
  
  # this loops over the list of models which have to be sorted
  # the loops starts right after the best model
  loop_sequence <- (best_pos + 1L):length( lda_models )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  for ( i_mod in loop_sequence ) {
    # for ( i_mod in seq_along( lda_models ) ) {
    i_pos <- i_mod - .optimal_model + 1L
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
    for ( j in 1L:.optimal_model ) {
      # find the max cosine similarity and return the value and its index
      maxval <- max( CosSim[ j, ] )
      maxsim <- which.max( CosSim[ j, ] )
      BestTopics[ (current_k - .optimal_model), j ] <- maxsim
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
  }
  
  Chi_K[ , 2L ] <- apply( Chi_k, 1L, sum ) / .optimal_model
  # not sure if we want to divide by 1:108 or by the real k given
  # by topics over optimal
  Chi_K[ , 3L ] <- Chi_K[ , 2L ] / Chi_K[ , 1L ]
  colnames( Chi_K ) <- c( "dfs", "chisq", "chisq_std" )
  
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  if ( convert ) {
    Chi_K <- as.data.frame( Chi_K )
    Chi_K$topic <- Chi_K$dfs + .optimal_model
    Chi_K[ , c( 4L, 1L:3L ) ]
  } else {
    Chi_K <- as.data.table( Chi_K )
    Chi_K[ , topic := dfs + .optimal_model ]
    setcolorder( Chi_K, c( 4L, 1L:3L ) )
  }
  if ( do_plot ) {
    cat( "---\n" )
    cat( "Plotting...\n" )
    p1 <- ggplot( Chi_K ) +
      geom_line( aes( x = topic, y = chisq_std ), size = 0.8, color = "royalblue" ) +
      xlab( "Topics" ) + ylab( expression( paste( "TWW - ", chi^2 ) ) )
    print( p1 )
  }
  
  toc <- proc.time()
  runtime <- toc - tic
  cat( "---\n" )
  cat( "Function took:", runtime[ 3L ], "sec.\n" )
  cat( "---\n" )
  return( Chi_K[] )
}
