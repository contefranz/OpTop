if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "id_word", "id_doc", "word_prop",
                             "word_sum", "check", "topics",
                             ".", "chisquare", "chisquare_mod",
                             "row_cut", "chi_sum", "word_prop_hat",
                             "word_prop_hat_cum", "pval" ) )
}
#' Find the optimal number of topics from a pool of LDA models
#'
#' Implements a fast chi-square like test to detect the number of topics
#' that best describes the corpus estimated via Latent Dirichelet Allocation.
#'
#' @param lda_models A list of ordered LDA models as computed by
#' \code{\link[topicmodels]{LDA}}. The LDA models must be in ascending order
#' according to the number of topics.
#' @param word_proportions A \code{\link[data.table]{data.table}} giving the 
#' word proportions in a corpus as computed by \code{\link[OpTop]{word_proportions}}.
#' @param q Set a cutoff for important words as the quantile of the expected
#' cumulative probability of word weights. Default to 0.80, meaning that the 
#' function reaches 80\% of the distribution mass and leaves out the remaining
#' 20\%.
#' @param alpha The confidence level of test acceptance. Default to 0.05. 
#' See 'Details'.
#' @param do_plot Plot the chi-square statistic as a function of the number of 
#' topics. Default to \code{TRUE}.
#' @param convert Target convertion format. This version of \code{OpTop} supports
#' \code{\link[base]{data.frame}} and \code{\link[tibble]{tibble}}.
#' Default to \code{NULL} which returns a \code{\link[data.table]{data.table}}.
#' @details The function implements a Pearson chi-square statistic that exploits
#' the assumption that the distribution of words is multinomial. The test studies
#' the stability of a K-topic model which fully characterizes the corpus if the
#' observed and estimated word vectors are statistically indistinct.
#' 
#' The parameter \code{alpha} controls the confidence of the chi-square test. The
#' optimal model is selected the first time the chi-square statistic reaches
#' a p-value equal to \code{alpha}. In the event that the chi-square statistic
#' fails to reach \code{alpha}, the minimum chi-square statistic
#' is selected. A higher \code{alpha} resolves in selecting a model with less 
#' topics. You can force the algorithm to find the minimum chi-square statistic
#' by setting \code{alpha} equal to zero.
#' @return A \code{data.table} containing the following columns:
#'
#' \item{\code{topic}}{An integer giving the number of topics.}
#' \item{\code{chisq_std}}{A numeric giving the standardized chi-square.}
#' \item{\code{pval}}{A numeric giving the p-value of the test.}
#' @examples
#' \dontrun{
#' # Compute word proportions from a corpus objects
#' word_proportions <- word_proportions( corpus = data_corpus_inaugural,
#'                                       remove_document = TRUE,
#'                                       language = "en",
#'                                       source = "snowball" )
#'
#' test1 <- optimal_topic( lda_models = lda_list,
#'                         word_proportions = word_proportions,
#'                         q = 0.80,
#'                         alpha = 0.05 )
#' }
#' @seealso \code{\link[topicmodels]{LDA}} \code{\link[data.table]{data.table}}
#' @references Lewis, C. and Grossetti, F. (2019 - forthcoming):\cr
#' A Statistical Approach for Optimal Topic Model Identification.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}
#' @author Craig M. Lewis \email{craig.lewis@@owen.vanderbilt.edu}
#' @import data.table ggplot2
#' @importFrom tibble as_tibble
#' @export

optimal_topic <- function( lda_models, word_proportions,
                           q = 0.80, alpha = 0.05, 
                           do_plot = TRUE, 
                           convert = NULL ) {
  
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
  if ( !is.numeric( q ) ) {
    stop( "q must be a numeric" )
  }
  if ( !is.numeric( alpha ) ) {
    stop( "alpha must be a numeric" )
  }
  if ( !is.null( convert ) && !is.character( convert ) ) {
    stop( "When not NULL, convert must be either a \"data.frame\" or a \"tibble\"" )
  }
  
  
  tic <- proc.time()
  # compute the size of vocabulary detected in each document as:
  size_corpus <- nrow( word_proportions )
  # size_vocabulary <- nrow( word_proportions[ , .N, by = id_word ] )
  size_vocabulary <- word_proportions[ , max( id_word ) ]
  n_docs <- word_proportions[ , max( id_doc ) ]
  
  # final output table
  # regstats <- data.table()
  regstats <- matrix( NA_real_, nrow = 0, ncol = 4 )
  Chi_K <- data.table()
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  # we enter in looping over each model (j)
  for ( i_mod in seq_along( lda_models ) ) {
    # getting the document word weights --> gamma
    # dww <- as.data.table( lda_models[[ i_mod ]]@gamma )
    dww <- lda_models[[ i_mod ]]@gamma
    current_k <- ncol( dww )
    cat( "---\n" )
    cat( "# # # Processing LDA with k =", current_k, "\n" )
    # getting the term word weights --> beta
    # tww <- as.data.table( t( exp( lda_models[[ i_mod ]]@beta ) ) )
    tww <- t( exp( lda_models[[ i_mod ]]@beta ) )
    # adding row position to both objects
    # dww[ , id_doc := .I ]
    dww <- cbind( dww, 1L:nrow(dww) )
    # tww[ , id_word := .I ]
    tww <- cbind( tww, 1L:nrow(tww) )
    
    # looping over each document (k) in each model (j)
    # this is the loop that needs to be parallelized
    cat( "--> Processing documents\n" )
    for ( j_doc in 1L:n_docs ) {
      # subsetting word proportions based on id_doc
      # prop <- word_proportions[ id_doc == j_doc ]
      prop <- word_proportions[ .(j_doc) ]
      # subsetting dww according to id_doc
      # dwwj_doc <- as.matrix( dww[ j_doc ] )
      dwwj_doc <- dww[ j_doc, ]
      
      # casting N x K matrix
      # dww_j_doc <- matrix( data = dwwj_doc,
      #                      ncol = ncol( dwwj_doc ),
      #                      nrow = size_vocabulary,
      #                      byrow = TRUE )
      dww_j_doc <- matrix( data = dwwj_doc,
                           ncol = length( dwwj_doc ),
                           nrow = size_vocabulary,
                           byrow = TRUE )
      
      # this avoids the use of j index which does not match with matlab code
      # in matlab j loops over k_start -> k_end
      # here starts from 1 up to the latest model
      sub_dww_j_doc <- dww_j_doc[ , 1L:( ncol(dww_j_doc) - 1L ) ]
      # sub_tww <- as.matrix( tww[ , 1L:( ncol(tww) - 1L ) ] )
      sub_tww <- tww[ , 1L:( ncol(tww) - 1L ) ]
      
      # dot product --> element-wise multiplication
      tww_dww <- sub_dww_j_doc * sub_tww
      
      # this returns a vector...maybe we want a matrix
      X <- base::rowSums( tww_dww )
      # BestPair <- data.table( prop[ , .( word_prop ) ],
      #                         word_prop_hat = X )
      BestPair <- cbind( prop[ , word_prop ],
                         X )
      # setorder( BestPair, -word_prop_hat )
      BestPair <- BestPair[ order(-BestPair[ , 2L ] ), ]
      # compute the cumlative probability over estimations
      # BestPair[ , word_prop_hat_cum := cumsum( word_prop_hat ) ]
      BestPair <- cbind( BestPair, cumsum( BestPair[ , 2L ] ) )
      n_BP <- nrow( BestPair )
      p_BP <- ncol( BestPair )
      # stop when you reach q
      # AggBestPair <- BestPair[ round( word_prop_hat_cum, 4L ) <= q ]
      AggBestPair <- BestPair[ which( round(BestPair[ , 3L ], 4L) <= q ), ]
      icut <- nrow( AggBestPair )
      # lowest_estimates <- BestPair[ (icut + 1L):n_BP, lapply( .SD, sum ) ]
      lowest_estimates <- apply( BestPair[ (icut + 1L):n_BP, ], 2L, sum )
      # AggBestPair <- rbindlist( list( AggBestPair, lowest_estimates ) )
      AggBestPair <- rbind( AggBestPair, lowest_estimates )
      # numerator <- ( AggBestPair[ , word_prop ] - AggBestPair[ , word_prop_hat ] )^2L
      # denominator <- AggBestPair[ , word_prop_hat ]
      numerator <- ( AggBestPair[ , 1L ] - AggBestPair[ , 2L ] )^2L
      denominator <- AggBestPair[ , 2L ]
      chi_sq_fit <- icut * sum( numerator / denominator )
      
      # column chisquare_mod is just a placeholder here
      # this is to avoid the duplication of regstats in the outer loop
      # regout <- data.table( topics = current_k,
      #                       id_doc = j_doc,
      #                       chisquare = chi_sq_fit,
      #                       row_cut = icut )
      regout <- cbind( current_k, j_doc, chi_sq_fit, icut )
      regstats <- rbind( regstats, regout )
      
    }
    
    # chi_out <- regstats[ topics == current_k ]
    chi_out <- regstats[ which( regstats[ , 1L ] == current_k ) , ]
    # sum_i_mod <- chi_out[ , .( chi_sum = sum( chisquare ), cut = sum( row_cut ) ) ]
    sum_i_mod <- cbind( sum( chi_out[ , 3L ] ), sum( chi_out[ , 4L ] ) )
    # temp <- data.table( topic = current_k,
    #                     chisq_std = sum_i_mod[ 1L, chi_sum ] / sum_i_mod[ 1L, cut ],
    #                     pval = NA_real_
    # )
    temp <- cbind( current_k, sum_i_mod[ , 1L ] / sum_i_mod[ , 2L ] )
    
    # temp[ , pval := pchisq( chisq_std, df = 1L ) ]
    temp <- cbind( temp, pchisq( temp[ , 2L ], df = 1L ) )
    Chi_K <- rbind( Chi_K, temp )
    
  }
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  cat( "---\n" )
  
  setnames( Chi_K, old = names( Chi_K ), c( "topic", "chisq_std", "pval" ) )  
  
  global_min <- Chi_K[ , .SD[ which.min( chisq_std ) ] ]
  alpha_min <- Chi_K[ pval <= alpha ][ 1L ]
  if ( alpha == 0 || all( is.na( alpha_min ) ) ) {
    cat( "Optimal model found by global minimum\n" )
    cat( "Optimal model has", global_min$topic, "topics\n" )
    best_topic <- global_min
  } else {
    if ( global_min$topic > alpha_min$topic ) {
      cat( "Optimal model found by significance level of", alpha, "\n" )
      cat( "Optimal model has", alpha_min$topic, "topics\n" )
      best_topic = alpha_min
    } else {
      cat( "Optimal model found by global minimum\n" )
      cat( "Optimal model has", global_min$topic, "topics\n" )
      best_topic <- global_min
    }
  }
  
  if ( do_plot ) {
    cat( "Plotting...\n" )
    x_min <- best_topic$topic
    y_min <- best_topic$chisq_std
    p1 <- ggplot( Chi_K ) +
      geom_line( aes( x = topic, y = chisq_std ), size = 0.8, color = "royalblue" ) +
      geom_hline( yintercept = y_min, color = "black", linetype = 2L ) +
      geom_vline( xintercept = x_min, color = "black", linetype = 2L ) +
      geom_point( aes( x = x_min, y_min ), color = "red", shape = 4L, size = 4L ) +
      scale_y_continuous( breaks = seq( 0, max( Chi_K$chisq_std ), by = 0.5 ) ) +
      xlab( "Topics" ) + ylab( expression( bold( chi^2 ) ) ) +
      ggtitle( "Optimal Topic Plot" ) +
      theme_OpTop
    print( p1 )
  }
  
  if ( !is.null( convert ) ) {
    cat( "Converting to", convert, "\n" )
    if ( convert == "data.frame" ) {
      setDF( Chi_K )
    } else if ( convert == "tibble" ) {
      Chi_K <- as_tibble( Chi_K )
    }
  }
  
  toc <- proc.time()
  runtime <- toc - tic
  cat( "---\n" )
  cat( "Function took:", runtime[ 3L ], "sec.\n" )
  cat( "---\n" )
  return( Chi_K )
}
