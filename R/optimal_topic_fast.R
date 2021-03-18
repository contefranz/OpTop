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
#' It is recommended that \code{weighted_dfm} has element names consistent with the ones detected
#' by \code{\link[topicmodels]{LDA}}. See 'Details'.
#' @param q Set a cutoff for important words as the quantile of the expected
#' cumulative probability of word weights. Default to 0.80, meaning that the 
#' function reaches 80\% of the distribution mass and leaves out the remaining
#' 20\%.
#' @param alpha The confidence level of test acceptance. Default to 0.05. 
#' See 'Details'.
#' @param ncores Control the number of cores to use through \code{\link[parallel]{detectCores}}.
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
#' has specific element names indicating the original names of the documents as defined
#' in the \code{\link[quanteda]{corpus}}. These element names can be extracted with 
#' \code{\link[quanteda]{docid}}\code{(weighted_dfm)}. 
#' If, for any reason, the function \code{\link[topicmodels]{LDA}} fails 
#' to estimate the requested \code{k} topics over certain documents, then \code{optimal_topic} 
#' takes care of that by ensuring that there is a perfect match between the documents found in 
#' \code{weighted_dfm} and the ones contained in \code{lda_models}. 
#' If \code{weighted_dfm} does not contain any meaningful name to be matched with \code{lda_models}, 
#' for instance if the whole vector is full of \code{FALSE}, then \code{optimal_topic} stops 
#' with an error because most likely there is something wrong. If \code{optimal_topic} finds
#' few documents that are not present in \code{lda_models}, then it removes them from the input
#' \code{weighted_dfm} in order to achieve a perfect match. 
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
#' @importFrom parallel detectCores clusterEvalQ clusterExport parLapply
#' @export

optimal_topic_fast <- function( lda_models, weighted_dfm,
                                q = 0.80, alpha = 0.05, 
                                ncores = NULL,
                                do_plot = TRUE, 
                                convert = NULL ) {
  
  if ( !is.list( lda_models ) ) {
    stop( "lda_models must be a list" )
  }
  if ( length( lda_models ) == 1L ) {
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
  if ( is.null( ncores ) ) {
    ncores = detectCores()
  } else {
    if ( !is.numeric( ncores ) ) {
      stop("ncores must be a positive integer")
    }
    clust = makeCluster(ncores)
  }
  
  tic <- proc.time()
  # compute the number of docs and features in the vocabulary
  docs <- as.character( docid( weighted_dfm ) )
  n_docs <- length( docs )
  n_features <- nfeat( weighted_dfm )
  features = featnames( weighted_dfm )
  
  # final output table
  regstats <- matrix( NA_real_, nrow = 0, ncol = 4 )
  Chi_K <- data.table()
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  cat("Using", ncores, "cores\n" )
  clusterEvalQ(clust, c(library(data.table), library(quanteda)))
  clusterExport( clust, 
                 c("lda_models", "weighted_dfm", "docs", "n_docs", "n_features", "features"), 
                 envir = environment())
  
  chisqlist = parLapply(cl = clust, X = seq_along( lda_models ), fun = function( i_mod ) {
    # we enter in looping over each model (j)
    # getting the document word weights --> gamma
    dww <- lda_models[[ i_mod ]]@gamma
    current_k <- ncol( dww )
    
    doc_check <- docs %in% lda_models[[ i_mod ]]@documents
    
    if ( all(doc_check) ) {
      cat("Found perfect match between LDA-documents and weighted_dfm\n" )
    } else {
      id_toremove <- which( doc_check == FALSE )
      if ( length( id_toremove ) < length( doc_check ) ) {
        cat("Removing unmatched documents\n" )
        toremove <- docs[ id_toremove ]
        weighted_dfm <- weighted_dfm[ -id_toremove, ]
      } else {
        stop("Document matching went really wrong. Check docs in both weighted_dfm and in LDA@documents")
      }
    }
    
    # if ( !all(doc_check) ) {
    #   stop("All documents are missing in lda_models. Something is wrong...")
    # } else if ( all(doc_check) ) {
    #   cat("Found perfect match between LDA-documents and weighted_dfm\n" )
    # } else {
    #   cat("Removing unmatched documents\n" )
    #   id_toremove <- which( doc_check == FALSE )
    #   if ( length( id_toremove ) < 1L ) {
    #     stop("Document matching went really wrong. Check docs in both weighted_dfm and in LDA@documents")
    #   }
    #   toremove <- docs[ id_toremove ]
    #   weighted_dfm <- weighted_dfm[ -id_toremove, ]
    # }
    # 
    # getting the term word weights --> beta
    tww <- t( exp( lda_models[[ i_mod ]]@beta ) )
    # adding row position to both objects
    
    tww_dww_long = apply(dww, 1L, function(x) 
      matrix(x, nrow = n_features, ncol = current_k, byrow = TRUE) * tww) #done yay!!!!
    
    tww_dww = lapply( seq_len( ncol( tww_dww_long ) ), function( x ) 
      matrix( tww_dww_long[ , x ], ncol = current_k ) )
    
    X = lapply( tww_dww, function( x ) {
      out = data.table( feature = features, 
                        X = base::rowSums(x) ) 
      setorder( out, -X )
      out[ , cumX := cumsum( X ) ]
      return( out )
    }
    )
    BestCut = lapply( X, function( x ) {
      highest = x[ cumX <= q ]
      lowest = x[ cumX > q, .( feature = NA_character_, X = sum( X ), cumX = sum( cumX ) ) ]
      out = rbindlist( list( highest, lowest ) )
    } )
    
    chi_sq_fit = lapply( seq_along(BestCut), function( i_doc ) {
      
      keep_features = na.omit(BestCut[[ i_doc ]][ , feature ] ) 
      current_dfm = dfm_subset( weighted_dfm, doc_id == docs[ i_doc ] )
      current_dfm = dfm_keep( current_dfm, pattern = keep_features, valuetype = "fixed" )
      current_features = featnames( current_dfm )
      dt_prop = convert( current_dfm, to = "data.frame")
      setDT(dt_prop)
      dt_prop[ , doc_id := NULL ]  
      dt_prop = data.table::transpose(dt_prop)
      dt_prop[ , feature := current_features ]
      res_prop = 1 - sum(dt_prop$V1)
      
      setkey( dt_prop, feature )
      setkey( BestCut[[ i_doc ]], feature )
      BestCut[[ i_doc ]] = dt_prop[ BestCut[[ i_doc ]] ]
      BestCut[[ i_doc ]][ is.na(V1), V1 := res_prop ]
      
      numerator = BestCut[[ i_doc ]][ , (V1 - X)^2 ]
      denominator = BestCut[[ i_doc ]][ , X ]
      icut = length(current_features)
      chi_sq_fit <- data.table( topic = current_k, 
                                doc_id = i_doc, 
                                chi_sq_fit = icut * sum( numerator / denominator ),
                                icut = icut
      )
    } 
    )
    
  }
  )
  stopCluster(clust)
  
  regstats = rbindlist(lapply( chisqlist, rbindlist))
  Chi_K = regstats[ , .( sum(chi_sq_fit), sum(icut) ), by = topic]
  Chi_K[ , OpTop := V1/V2 ]  
  Chi_K[ , pval := pchisq( OpTop, df = 1L ) ]
  Chi_K = Chi_K[ , .( topic, OpTop, pval ) ]
  
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  cat( "---\n" )

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
  # return( Chi_K )
}
