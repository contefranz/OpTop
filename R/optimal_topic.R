if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "id_word", "id_doc", "weighted_dfm", 
                             "word_prop", "document",
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
#' @param weighted_dfm A weighted \code{\link[quanteda]{dfm}} containing word proportions.
#' It is recommended that \code{weighted_dfm} has the \code{\link[quanteda]{docvar}} 
#' "doc_id" with original document names. See 'Details'.
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
#' To ensure a complete matching between the set of LDA models specified 
#' through \code{lda_models}, we strongly recommend that the corresponding \code{weighted_dfm} 
#' has a specific \code{docvar} indicating the original names of the documents as defined
#' in the \code{corpus}. If, for any reason, the function \code{\link[topicmodels]{LDA}} fails 
#' to estimate the requested \code{k} topics, then \code{optimal_topic} takes care of that 
#' by ensuring that there is a perfect overlap between the documents found in \code{weighted_dfm}
#' and the ones contained in \code{lda_models}. If \code{weighted_dfm} does not contain such 
#' \code{docvar}, then \code{optimal_topic} throws a warning and assumes a complete match between
#' the documents and those estimated by \code{\link[topicmodels]{LDA}}.
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
#'
#' \item{\code{topic}}{An integer giving the number of topics.}
#' \item{\code{OpTop}}{A numeric giving the standardized chi-square.}
#' \item{\code{pval}}{A numeric giving the p-value of the test.}
#' @examples
#' \dontrun{
#' # Compute word proportions from a corpus objects
#' weighted_dfm <- weighted_dfm( corpus = data_corpus_inaugural,
#'                                       remove_document = TRUE,
#'                                       language = "en",
#'                                       source = "snowball" )
#'
#' test1 <- optimal_topic( lda_models = lda_list,
#'                         weighted_dfm = weighted_dfm,
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
#' @importFrom quanteda ndoc nfeat is.dfm
#' @export

optimal_topic <- function( lda_models, weighted_dfm,
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
  if( !is.dfm( weighted_dfm ) ) {
    stop( "weighted_dfm must be a dfm" )
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
  # compute the number of docs and features in the vocabulary
  n_docs <- ndoc( weighted_dfm )
  n_features <- nfeat( weighted_dfm )
  
  # final output table
  regstats <- matrix( NA_real_, nrow = 0, ncol = 4 )
  Chi_K <- data.table()
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  # we enter in looping over each model (j)
  for ( i_mod in seq_along( lda_models ) ) {
    
    # getting the document word weights --> gamma
    dww <- lda_models[[ i_mod ]]@gamma
    current_k <- ncol( dww )
    cat( "---\n" )
    cat( "# # # Processing LDA with k =", current_k, "\n" )
    
    # if ( "doc_id" %in% names( weighted_dfm ) ) {
    #   cat( "Checking which documents have been estimated by LDA\n" )
    #   docs <- unique( weighted_dfm$doc_id )
    #   doc_check <- docs %in%
    #     lda_models[[ i_mod ]]@documents
    #   setkey( weighted_dfm, doc_id )
    #   weighted_dfm <- weighted_dfm[ !.( docs[ !doc_check ] ) ]
    #   n_docs <- uniqueN( weighted_dfm$id_doc )
    #   if ( length( which( doc_check == FALSE ) ) > 0 ) {
    #     weighted_dfm[ , id_doc := .GRP, by = doc_id ]
    #   }
    # } else {
    #   warning( "doc_id is not in weighted_dfm. ",
    #   "Assuming that there is a complete overlap between documents in ",
    #   "the corpus/dfm and those estimated by LDA()." )
    #   n_docs <- uniqueN( weighted_dfm$id_doc )
    # }
    
    # getting the term word weights --> beta
    tww <- t( exp( lda_models[[ i_mod ]]@beta ) )
    # adding row position to both objects
    dww <- cbind( dww, 1L:nrow(dww) )
    tww <- cbind( tww, 1L:nrow(tww) )
    
    # looping over each document (k) in each model (j)
    # this is the loop that needs to be parallelized
    cat( "--> Processing documents\n" )
    for ( j_doc in 1L:n_docs ) {
      # subsetting word proportions based on id_doc
      prop <- matrix(weighted_dfm[ j_doc, ])
      # subsetting dww according to id_doc
      dwwj_doc <- dww[ j_doc, ]
      
      # casting N x K matrix
      dww_j_doc <- matrix( data = dwwj_doc,
                           ncol = length( dwwj_doc ),
                           nrow = n_features,
                           byrow = TRUE )
      
      # this avoids the use of j index which does not match with matlab code
      # in matlab j loops over k_start -> k_end
      # here starts from 1 up to the latest model
      sub_dww_j_doc <- dww_j_doc[ , 1L:( ncol(dww_j_doc) - 1L ) ]
      sub_tww <- tww[ , 1L:( ncol(tww) - 1L ) ]
      
      # dot product --> element-wise multiplication
      tww_dww <- sub_dww_j_doc * sub_tww
      
      # this returns a vector...maybe we want a matrix
      X <- base::rowSums( tww_dww )
      BestPair <- cbind( prop, X )
      BestPair <- BestPair[ order(-BestPair[ , 2L ] ), ]
      # compute the cumlative probability over estimations
      BestPair <- cbind( BestPair, cumsum( BestPair[ , 2L ] ) )
      n_BP <- nrow( BestPair )
      p_BP <- ncol( BestPair )
      # stop when you reach q
      AggBestPair <- BestPair[ which( round(BestPair[ , 3L ], 4L) <= q ), ]
      icut <- nrow( AggBestPair )
      lowest_estimates <- apply( BestPair[ (icut + 1L):n_BP, ], 2L, sum )
      AggBestPair <- rbind( AggBestPair, lowest_estimates )
      numerator <- ( AggBestPair[ , 1L ] - AggBestPair[ , 2L ] )^2L
      denominator <- AggBestPair[ , 2L ]
      chi_sq_fit <- icut * sum( numerator / denominator )
      
      # column chisquare_mod is just a placeholder here
      # this is to avoid the duplication of regstats in the outer loop
      regout <- cbind( current_k, j_doc, chi_sq_fit, icut )
      regstats <- rbind( regstats, regout )
      
    }
    
    chi_out <- regstats[ which( regstats[ , 1L ] == current_k ) , ]
    sum_i_mod <- cbind( sum( chi_out[ , 3L ] ), sum( chi_out[ , 4L ] ) )
    temp <- cbind( current_k, sum_i_mod[ , 1L ] / sum_i_mod[ , 2L ] )
    temp <- cbind( temp, pchisq( temp[ , 2L ], df = 1L ) )
    Chi_K <- rbind( Chi_K, temp )
    
  }
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  cat( "---\n" )
  
  setnames( Chi_K, old = names( Chi_K ), c( "topic", "OpTop", "pval" ) )  
  
  global_min <- Chi_K[ , .SD[ which.min( OpTop ) ] ]
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
    y_min <- best_topic$OpTop
    p1 <- ggplot( Chi_K ) +
      geom_line( aes( x = topic, y = OpTop ), size = 0.8, color = "royalblue" ) +
      geom_hline( yintercept = y_min, color = "black", linetype = 2L ) +
      geom_vline( xintercept = x_min, color = "black", linetype = 2L ) +
      geom_point( aes( x = x_min, y_min ), color = "red", shape = 4L, size = 4L ) +
      xlab( "Topics" ) + ylab( expression(OpTop[J]^{"K"}) ) +
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
