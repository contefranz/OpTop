if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "id_word", "id_doc", "word_prop",
                             "word_sum", "check", "topics",
                             ".", "chisquare", "chisquare_mod",
                             "row_cut", "chi_sum" ) )
}
#' Find the optimal number of topics from a pool of LDA models
#'
#' Implements a fast chi-square like test to detect the number of topics
#' that best describes the corpus.
#'
#' @param lda_models A list of ordered LDA models as computed by
#' \code{\link[topicmodels]{LDA}}. The LDA models must be in ascending order
#' according to the number of topics.
#' @param word_proportions A \code{data.table} giving the word proportions in a corpus
#' as computed by \code{\link[OpTop]{word_proportions}}.
#' @param threshold Set a cutoff for the relative importance of words.
#' Default to 0.00075.
#' @param alpha Probability at which compute the quantiles of the chi-square test.
#' Default to 0.01.
#' @param q_type Select the quantile algorithm as in \code{\link[stats]{quantile}}.
#' Default to 5 for consistency with Matlab. In later releases, it will be
#' replaced with 7, which is the \code{R} default.
#' @param convert Convert the output to a \code{data.frame}.
#' Default to \code{FALSE}.
#' @details The function implements a Pearson chi-square statistic that exploits
#' the assumption that the distribution of words is multinomial. The test studies
#' the stability of a K-topic model which fully characterizes the corpus if the
#' observed and estimated word vectors are statistically indistinct.
#' @return A \code{data.table} containing the following columns:
#'
#' \item{\code{topic}}{An integer giving the number of topics.}
#' \item{\code{chi_sum}}{A numeric giving the overall sum of the chi-square.}
#' \item{\code{cut}}{An integer giving the row position of minimum word
#' importance.}
#' \item{\code{chi_std}}{A numeric giving the standardized chi-square.}
#' \item{\code{chi_stdw}}{A numeric giving the standardized and weighted chi-square.}
#' @examples
#' \dontrun{
#' # Compute word proportions from a corpus objects
#' word_proportions <- word_proportions( corpus = data_corpus_inaugural,
#'                                      remove_document = TRUE,
#'                                      language = "en",
#'                                      source = "snowball" )
#'
#' test1 <- optimal_topic( lda_models = lda_list,
#'                        word_proportions = word_proportions,
#'                        threshold = 0.00075,
#'                        alpha = 0.01,
#'                        q_type = 5 )
#' }
#' @seealso \code{\link[topicmodels]{LDA}} \code{\link[data.table]{data.table}}
#' @references Lewis, C. and Grossetti, F. (2019 - forthcoming):\cr
#' A Statistical Approach for Optimal Topic Model Identification.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}.
#' @import data.table
#' @importFrom stats quantile
#' @export

optimal_topic = function( lda_models, word_proportions,
                          threshold = 0.00075, alpha = 0.01,
                          q_type = 5L, convert = FALSE ) {

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
  if ( !is.numeric( alpha ) ) {
    stop( "alpha must be a numeric" )
  }
  if ( !is.numeric( q_type ) || q_type < 1L || q_type > 9L ) {
    stop( "q_type must be an integer between 1 and 9" )
  }
  if ( !is.logical( convert ) ) {
    stop( "convert must be either TRUE or FALSE" )
  }

  tic <- proc.time()
  # compute the size of vocabulary detected in each document as:
  size_corpus <- nrow( word_proportions )
  size_vocabulary <- nrow( word_proportions[ , .N, by = id_word ] )
  n_docs <- nrow( word_proportions[ , .N, by = id_doc ] )

  # final output table
  regstats <- data.table()
  Chi_K <- data.table()
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  # we enter in looping over each model (j)
  for ( i_mod in seq_along( lda_models ) ) {
    # getting the document word weights --> gamma
    dww <- as.data.table( lda_models[[ i_mod ]]@gamma )
    current_k <- ncol( dww )
    cat( "---\n" )
    cat( "# # # Processing LDA with k =", current_k, "\n" )
    # getting the term word weights --> beta
    tww <- as.data.table( t( exp( lda_models[[ i_mod ]]@beta ) ) )
    # adding row position to both objects
    dww[ , id_doc := .I ]
    tww[ , id_word := .I ]

    # looping over each document (k) in each model (j)
    # this is the loop that needs to be parallelized
    cat( "--> Processing documents\n" )
    for ( j_doc in 1L:n_docs ) {
      # subsetting word proportions based on id_doc
      prop <- word_proportions[ id_doc == j_doc ]
      # subsetting dww according to id_doc
      dwwj_doc <- as.matrix( dww[ j_doc ] )

      # casting N x K matrix
      dww_j_doc <- matrix( data = dwwj_doc,
                      ncol = ncol( dwwj_doc ),
                      nrow = size_vocabulary,
                      byrow = TRUE )

      # this avoids the use of j index which does not match with matlab code
      # in matlab j loops over k_start -> k_end
      # here starts from 1 up to the latest model
      sub_dww_j_doc <- dww_j_doc[ , 1L:( ncol(dww_j_doc) - 1L ) ]
      sub_tww <- as.matrix( tww[ , 1L:( ncol(tww) - 1L ) ] )

      # dot product --> element-wise multiplication
      tww_dww <- sub_dww_j_doc * sub_tww

      # this returns a vector...maybe we want a matrix
      X <- base::rowSums( tww_dww )
      BestPair <- data.table( prop[ , .( word_prop ) ], word_sum = X )
      setorder( BestPair, -word_sum )
      # take the minimum of estimated value
      # save the minimum value and the location of the minimum (i.e. row)
      BestPair[ , check := abs( word_sum - threshold ) ]
      pct <- BestPair[ which.min( check ), check ]
      icut <- which.min( BestPair[ , check ] )
      BestPair[ , check := NULL ]
      n_BP <- nrow( BestPair )
      p_BP <- ncol( BestPair )
      AggBestPair <- BestPair[ 1L:icut ]
      lowest_estimates <- BestPair[ (icut + 1L):n_BP, lapply( .SD, sum ) ]
      AggBestPair <- rbindlist( list( AggBestPair, lowest_estimates ) )
      numerator <- ( AggBestPair[ , word_prop ] - AggBestPair[ , word_sum ] )^2L
      denominator <- AggBestPair[ , word_sum ]
      chi_sq_fit <- icut * sum( numerator / denominator )

      # column chisquare_mod is just a placeholder here
      # this is to avoid the duplication of regstats in the outer loop
      regout <- data.table( topics = current_k,
                           id_doc = j_doc,
                           chisquare = chi_sq_fit,
                           row_cut = icut )
      regstats <- rbind( regstats, regout )

    }

    chi_out <- regstats[ topics == current_k ]
    # NOTE: Matlab uses quantile type = 5 ?!?!
    # R by default uses type = 7
    data_min <- stats::quantile( chi_out[ , chisquare ],
                                probs = alpha,
                                type = q_type )
    data_max <- stats::quantile( chi_out[ , chisquare ],
                                probs = 1L - alpha,
                                type = q_type )

    # update chisquare with quantiles but keep trace of the original one
    chi_out[ , chisquare_mod := chisquare ]
    chi_out[ chisquare < data_min, chisquare_mod := data_min ]
    chi_out[ chisquare > data_max, chisquare_mod := data_max  ]
    sum_i_mod <- chi_out[ , .( chi_sum = sum( chisquare ), cut = sum( row_cut ) ) ]
    temp <- data.table( topic = current_k,
                       sum_i_mod,
                       chi_std = sum_i_mod[ 1L, chi_sum ] / sum_i_mod[ 1L, cut ],
                       chi_stdw = chi_out[ , sum( chisquare_mod ) ] / sum_i_mod[ 1L, cut ] )
    Chi_K <- rbind( Chi_K, temp )

  }
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  if ( convert ) {
    cat( "Converting to data.frame\n" )
    setDF( Chi_K )
  }

  toc <- proc.time()
  runtime <- toc - tic
  cat( "---\n" )
  cat( "Function took:", runtime[ 3L ], "sec.\n" )
  cat( "---\n" )
  return( Chi_K )
}
