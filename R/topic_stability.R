if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "chisq", "chi_std", "topic", "pchisq", "dfs", "chisq_std" ) )
}
#' Compute topic stability for over-optimal topic specifications
#'
#' Implements fast chi-square like test to evaluate the stability of redundant
#' topics.
#' 
#' @inheritParams optimal_topic
#' @param optimal_model A number corresponding to the optimal topic model.
#' @param alpha Alpha level to identify informative words from the Cumulative
#' Distribution Function over the cosine similarities in the Topic Word Weights
#' matrix. Default to 0.05.
#' @param do_plot Plot the chi-square statistic as a function of the number of 
#' topics. Default to \code{TRUE}. The horizontal dot-dashed line represents
#' the significance level according to \code{alpha}.
#' @details This function implements Test 3 as defined in 
#' Lewis and Grossetti (2019). Test 3 evaluates the aggregated 
#' stability of over-optimal topic specifications by summing each 
#' point-wise contribution. See 'Value' to understand how \code{topic_stability} 
#' returns the results.
#' @return A `data.table` containing the following columns:
#'
#' \item{\code{topic}}{An integer giving the number of topics.}
#' \item{\code{df}}{An integer giving the degrees of freedom.}
#' \item{\code{chisq}}{A numeric giving the chi-square statistic.}
#' @examples
#' \dontrun{
#' test2 <- topic_stability( lda_models = lda_list,
#'                           optimal_model = test1,
#'                           q = 0.00075,
#'                           alpha = 0.05 )
#' }
#' @seealso \code{\link[topicmodels]{LDA}} \code{\link[data.table]{data.table}}
#' @references Lewis, C. and Grossetti, F. (2019 - forthcoming):\cr
#' A Statistical Approach for Optimal Topic Model Identification.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}
#' @author Craig M. Lewis \email{craig.lewis@@owen.vanderbilt.edu}
#' @import data.table ggplot2
#' @importFrom stats pchisq
#' @importFrom tibble as_tibble
#' @export

topic_stability <- function( lda_models, optimal_model,
                             q = 0.80, alpha = 0.05,
                             do_plot = TRUE ) {
  
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
  if ( !is.numeric( q ) ) {
    stop( "q must be a numeric" )
  }
  if ( !is.numeric( alpha ) ) {
    stop( "alpha must be a numeric" )
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
  
  #############################
  #############################
  # FROM HERE ON WE MOVE TO C++
  #############################
  # find the element corresponding to the best topic
  if ( !is.LDA_VEM( optimal_model ) ) {
    # extracting information from best model
    tww_best <- t( exp( lda_models[[ best_pos ]]@beta ) )
  }
  # Normalizing by scaling by vector norms
  tww_best_norm <- .Call(`_OpTop_normalize_columns`, tww_best)
  l_models = length( lda_models )
  # pre-allocating required objects
  BestTopics    <- matrix( 0, nrow = (l_models - best_pos), ncol = optimal_model )
  Chi_k         <- matrix( 0, nrow = (l_models - best_pos), ncol = (optimal_model + 1L) )
  model_list    <- unname( sapply( lda_models, function(x) x@k ) )
  Chi_k[ , 1L ] <- model_list[ model_list > optimal_model ]
  Chi_K         <- matrix( 0, nrow = (l_models - best_pos), ncol = 3L )
  Chi_K[ , 1L ] <- 1L:nrow( Chi_K )
  p             <- matrix( 0, nrow = (l_models - best_pos), ncol = (optimal_model + 1L) )
  CosSimHold    <- NULL

  # this loops over the list of models which have to be sorted
  # the loops starts right after the best model
  loop_sequence <- (best_pos + 1L):length( lda_models )
  min_loop <- min( loop_sequence )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  for ( i_mod in loop_sequence ) {
    i_pos <- i_mod - min_loop + 1L
    # reading tww
    tww <- t( exp( lda_models[[ i_mod ]]@beta ) )
    current_k <- ncol( tww )
    n_tww <- nrow( tww )
    p_tww <- ncol( tww )
    cat( "---\n" )
    cat( "# # # Processing LDA with k =", current_k, "\n" )
    
    # Normalizing by scaling by vector norms
    tww_norm <- .Call(`_OpTop_normalize_columns`, tww)
    
    # this is the matrix multiplication to get the Cosine Similarities
    # between the best set of topics and the targe ones
    CosSim <- t( tww_best_norm ) %*% tww_norm
    # transpose the result so that N reflects the number of topics under
    # investigation
    minCosSim <- t( CosSim )
    CumDF <- 0
    
    # and here begins the second loop
    cat( "--> Finding topics with highest cosine similarity\n" )
    for ( j in 1L:optimal_model ) {
      # find the max cosine similarity and return the value and its index
      maxval <- max( CosSim[ j, ] )
      maxsim <- which.max( CosSim[ j, ] )
      BestTopics[ i_pos, j ] <- maxsim
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
      n_BP <- nrow( BestPair )
      p_BP <- ncol( BestPair )
      
      # stop when you reach q
      AggBestPair <- BestPair[ BestPair[ , 3L ] >= 1L - q, ]
      icut <- nrow( AggBestPair )
      lowest_estimates <- apply( BestPair[ (icut + 1L):n_BP, ], 2L, sum )
      
      AggBestPair <- rbind( BestPair[ 1L:icut, ], unname( lowest_estimates ) )
      numerator <- ( AggBestPair[ , 1L ] - AggBestPair[ , 2L ] )^2L
      denominator <- AggBestPair[ , 1L ]
      Chi_k[ i_pos, j + 1L ] <- icut * sum( numerator / denominator )
      p[ i_pos, j ] <- pchisq( Chi_k[ i_pos, j + 1L ], icut )
    }
  }
  
  Chi_K[ , 2L ] <- apply( Chi_k, 1L, sum ) / optimal_model
  # not sure if we want to divide by 1:108 or by the real k given
  # by topics over optimal
  Chi_K[ , 3L ] <- Chi_K[ , 2L ] / optimal_model
  
  Chi_K <- cbind( model_list[ model_list > optimal_model ],
                  Chi_K[ , 3L ], 
                  pchisq( Chi_K[ , 3L ], df = 1L ) )
  colnames( Chi_K ) <- c( "topic", "chisq", "pval" )
  Chi_K <- as.data.table( Chi_K )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  prop_H0 = nrow( Chi_K[ pval >= alpha ] ) / nrow( Chi_K )
  if ( prop_H0 == 1 ) {
    cat( "Null hypothesis is always accepted at a level of", alpha, "\n" )
  } else {
    reject = Chi_K[ pval < alpha ]
    cat( "Null hypothesis is rejected at a level of", alpha, 
         "for the following models\n" )
    cat( "Number of topics:", 
         paste0( reject[ , topic ], collapse = ", " ), "\n" )
    cat( "Overall, topic stability is achieved for ", 
         round( prop_H0 * 100, 2 ), "% of the models\n", sep = "" )
  }
  
  if ( do_plot ) {
    cat( "---\n" )
    cat( "Plotting...\n" )
    p1 <- ggplot( Chi_K ) +
      geom_hline( yintercept = qchisq( alpha, 1L ), linetype = 4 ) +
      geom_line( aes( x = topic, y = chisq ), size = 0.8, color = "royalblue" ) +
      xlab( "Topics" ) + ylab( expression( bold( chi^2 ) ) ) +
      ggtitle( "Topic Stability Plot" ) +
      theme_OpTop
    print( p1 )
  }
  
  toc <- proc.time()
  runtime <- toc - tic
  cat( "---\n" )
  cat( "Function took:", runtime[ 3L ], "sec.\n" )
  cat( "---\n" )
  return( Chi_K[] )
}
