if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "chi_sq_std", "df", ".", "pval_inform",
                             "chisq_inform_std", "pval_uninform", 
                             "chisq_uninform_std", "Fstat", "pval_Fstat",
                             "chisq_smooth" ) )
}
#' Compute aggregate document stability and F-test
#'
#' @description 
#' Detects informative and uninformative components to compute aggregate document
#' stability. Performs a chi-square test to evaluate document stability, 
#' Also, computes a F-test to further evaluate deviation from optimal model.
#' 
#' @inheritParams agg_topic_stability
#' @param weighted_dfm A weighted \code{\link[quanteda]{dfm}} containing word proportions.
#' It is recommended that \code{weighted_dfm} has the corresponding internal variable that can be
#' accessed with \code{docid}. See ?\code{\link[OpTop]{optimal_topic}} for more details.
#' @return A \code{data.table} containing the following columns:
#'
#' \item{\code{topic}}{An integer giving the number of topics.}
#' \item{\code{id_doc}}{An integer document id as given in the original corpus.}
#' \item{\code{chisq_inform_std}}{A numeric giving the standardized chi-square statistic
#' for the informative component.}
#' \item{\code{chisq_uninform_std}}{A numeric giving the standardized chi-square statistic
#' for the uninformative component.}
#' \item{\code{pval_inform}}{A numeric giving the p-value of the chi-square test
#' over the informative component.}
#' \item{\code{pval_uninform}}{A numeric giving the p-value of the chi-square test
#' over the uninformative component.}
#' \item{\code{Fstat}}{A numeric giving the standardized F statistic
#' of the ratio \code{chisq_inform_std}/\code{chisq_uninform_std}.}
#' \item{\code{pval_Fstat}}{A numeric giving the p-value of the F test.}
#' @examples
#'\dontrun{
#' test4 <- agg_document_stability( lda_models = lda_list,
#'                                  weighted_dfm = weighted_dfm,
#'                                  smoothed = TRUE, do_plot = TRUE )
#' }
#' @seealso \code{\link[topicmodels]{LDA}} \code{\link[data.table]{data.table}}
#' @references Lewis, C. and Grossetti, F. (2019 - forthcoming):\cr
#' A Statistical Approach for Optimal Topic Model Identification.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}.
#' @author Craig M. Lewis \email{craig.lewis@@owen.vanderbilt.edu}
#' @import data.table
#' @import grid
#' @importFrom quanteda ndoc nfeat is.dfm
#' @importFrom gridExtra marrangeGrob
#' @importFrom stats pf qf qchisq pchisq loess fitted.values
#' @export

agg_document_stability <- function( lda_models, weighted_dfm, 
                                    optimal_model, 
                                    q = 0.80, alpha = 0.05, 
                                    smoothed = TRUE, do_plot = TRUE ) {
  
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
  if ( !is.numeric( optimal_model ) || optimal_model < 2 ) {
    stop("optimal_model must be a number greater than 2")
  }
  if ( optimal_model == lda_models[[ length(lda_models ) ]]@k ) {
    message("Optimal model is already the last one in lda_models. There is nothing to compute above that.")
    return( NULL )
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
  if ( !is.logical( smoothed ) ) {
    stop( "smoothed must be either TRUE or FALSE" )
  }
  if ( .optimal_model == lda_models[[ length(lda_models ) ]]@k ) {
    message("Optimal model is already the last one in lda_models. There is nothing to compute above that.")
    return( NULL )
  }
  ###############################################################
  # THIS IS OLD CODE TO BE REMOVED ONCE THE CONVERSION IS STABLE
  ###############################################################
  # 
  # if ( ( !is.data.frame( optimal_model ) || !is.data.table( optimal_model ) ) &&
  #      !is.numeric( optimal_model ) && !is.LDA_VEM( optimal_model ) ) {
  #   stop( paste( "optimal_model must be either 1. an integer identifying",
  #                "the number of topics which best fits the corpus",
  #                "2. a data.table/data.frame as returned by .optimal_model()",
  #                "3. an LDA_VEM obect as computed by topicmodels::LDA()" ) )
  # }
  # if ( is.numeric( optimal_model ) ) {
  #   .optimal_model <- optimal_model
  # } else if ( is.data.table( optimal_model ) || is.data.frame( optimal_model ) ) {
  #   cat( "optimal_model is a data.table or a data.frame.",
  #        "Extracting information about optimal model...\n" )
  #   .optimal_model <- optimal_model[ which.min( OpTop ), topic ]
  # } else if ( is.LDA_VEM( optimal_model ) ) {
  #   cat( "optimal_model is a LDA_VEM object.", 
  #        "Extracting information about the optimal model...\n" )
  #   dtw_best <- optimal_model@gamma
  #   tww_best <- t( exp( optimal_model@beta ) )
  #   .optimal_model <- ncol( dtw_best )
  # }
  # cat( "best model has", .optimal_model, "topics\n" )
  ###############################################################
  
  tic <- proc.time()
  best_pos <- which( sapply( lda_models, function( x ) x@k ) == optimal_model )
  if ( length( best_pos ) == 0 ) {
    stop("optimal_model does not correspond to any topic number in lda_models")
  }
  
  ###############################################################
  # THIS IS OLD CODE TO BE REMOVED ONCE THE CONVERSION IS STABLE
  ###############################################################
  # compute the number of docs and features in the vocabulary
  # n_docs <- ndoc( weighted_dfm )
  # n_features <- nfeat( weighted_dfm )
  # 
  # k_end <- max( sapply( lda_models, function( x ) x@k ) )
  # best_pos <- which( sapply( lda_models, function( x ) x@k ) == .optimal_model )
  # if ( length( best_pos ) == 0 ) {
  #   stop( paste( "There is no optimal model in lda_models.",
  #                "This could be either due to a wrong specification of",
  #                "argument optimal_model or",
  #                "if optimal_model is a data.table, the optimal model cannot be found",
  #                "in the list lda_models." ) )
  # }
  ###############################################################
  
  # THE CODE ABOVE HAS BEEN REPLACED BY THE FOLLOWING
  # THIS IS WHAT WE DO IN OPTIMAL_TOPIC()
  docs = as.character( docid( weighted_dfm ) )
  n_docs = ndoc( weighted_dfm )
  n_features = nfeat( weighted_dfm )
  
  # get the list of documents to work on, by removing those which are not in the LDA models
  # ASSUMPTION: we assume that whenever the LDA fails to estimate topics in a given document, 
  # that document is dropped unconditionally on LDA specifications. That is, if we set k = 2 or
  # k = 10, the same document wil be dropped. Hence, the original loop over "lda_models" does not
  # make sense anymore. 
  # 
  # SOLUTION: we only check once and for all on the optimal topic model
  doc_check = docs %in% lda_models[[ best_pos ]]@documents
  if ( !all(doc_check) ) {
    id_toremove = which( doc_check == FALSE )
    if ( length( id_toremove ) < length( doc_check ) ) {
      cat("Removing unmatched documents\n" )
      weighted_dfm = weighted_dfm[ -id_toremove, ]
      # update number of docs
      n_docs = ndoc( weighted_dfm )
    } else {
      stop("Document matching went really wrong. Check docs in both weighted_dfm and in LDA@documents")
    }
  }  
  ##########################
  # C++ BEGINS HERE ! ! !
  ##########################
  
  # TO DO: this also needs to be converted in C++
  if ( !is.LDA_VEM( optimal_model ) ) {
    # extracting information from best model
    dtw_best <- lda_models[[ best_pos ]]@gamma
    tww_best <- t( exp( lda_models[[ best_pos ]]@beta ) )
  }
  K_fitval <- matrix( data = 0, nrow = n_features, ncol = n_docs )
  cat( "---\n" )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  for ( j_doc in 1L:n_docs ) {
    
    prop <- matrix(weighted_dfm[ j_doc, ])
    # subsetting dtw according to id_doc
    dtwj_doc <- dtw_best[ j_doc, ]
    # casting N x K matrix
    dtw_j_doc <- matrix( data = dtwj_doc,
                         ncol = .optimal_model,
                         nrow = n_features,
                         byrow = TRUE )
    fitval <- (dtw_j_doc * tww_best) %*% matrix( data = 1L, 
                                                 nrow = .optimal_model,
                                                 ncol = 1L )
    K_fitval[ , j_doc ] <- fitval
  }
  
  loop_sequence <- (best_pos + 1L):length( lda_models )  
  min_loop <- min( loop_sequence )
  k = 1L
  Chi_K <- data.table()
  for ( i_mod in loop_sequence ) {
    i_pos <- i_mod - min_loop + 1L
    dtw <- lda_models[[ i_mod ]]@gamma
    tww <- t( exp( lda_models[[ i_mod ]]@beta ) )
    current_k <- ncol( tww )
    n_tww <- nrow( tww )
    p_tww <- ncol( tww )
    cat( "---\n" )
    cat( "# # # Processing LDA with k =", current_k, "\n" )
    # fix the index i down below
    i <- i_pos
    out_in <- data.table()
    cat( "--> Processing documents\n" )
    for ( j_doc in 1L:n_docs ) {
      
      # subsetting dtw according to id_doc
      dtwj_doc <- dtw[ j_doc, ]
      # casting N x K matrix
      dtw_j_doc <- matrix( data = dtwj_doc,
                           ncol = p_tww,
                           nrow = n_features,
                           byrow = TRUE )
      inform <- dtw_j_doc[ , 1L:.optimal_model ] * tww[ , 1L:.optimal_model ]
      inform_norm <- ( inform %*% matrix( data = 1L, 
                                          nrow = .optimal_model, 
                                          ncol = 1L ) ) / 
        sum( inform %*% matrix( data = 1L, 
                                nrow = .optimal_model, 
                                ncol = 1L ) )
      
      uninform <- dtw_j_doc[ , (.optimal_model+1L):p_tww ] * 
        tww[ , (.optimal_model+1L):p_tww ]
      uninform_norm <- ( uninform %*% matrix( data = 1L, 
                                              nrow = current_k - .optimal_model, 
                                              ncol = 1L ) ) / 
        sum( uninform %*% matrix( data = 1L, 
                                  nrow = current_k - .optimal_model, 
                                  ncol = 1L ) )
      X <- cbind( K_fitval[ , j_doc ], inform_norm, uninform_norm )
      BestPair <- apply( X, 2L, function( x ) sort( x, decreasing = TRUE ) )
      BestPair <- cbind( BestPair, cumsum( BestPair[ , 1L ] ) )
      n_BP <- nrow( BestPair )
      p_BP <- ncol( BestPair )
      # stop when you reach q
      AggBestPair <- BestPair[ BestPair[ , 4L ] <= q, ]
      icut <- nrow( AggBestPair )
      lowest_estimates <- apply( BestPair[ (icut + 1L):n_BP, ], 2L, sum )
      AggBestPair <- rbind( BestPair[ 1L:icut, ], unname( lowest_estimates ) )
      
      chisq_inform <- (icut  + 1L) * sum( ( AggBestPair[ , 2L ] - AggBestPair[ , 1L ] )^2L /
                                            AggBestPair[ , 1L ] )
      chisq_uninform <- (icut + 1L) * sum( ( AggBestPair[ , 3L ] - AggBestPair[ , 1L ] )^2L /
                                             AggBestPair[ , 1L ] )
      out_doc <- data.table( topic = current_k,
                             id_doc = j_doc,
                             df = icut,
                             chisq_inform = chisq_inform,
                             chisq_uninform = chisq_uninform )
      out_in <- rbindlist( list( out_in, out_doc ) )
      
    }
    Chi_K <- rbindlist( list( Chi_K, out_in ) )
  }
  Chi_K[ , `:=` ( chisq_inform_std = chisq_inform / df,
                  chisq_uninform_std = chisq_uninform / df) ]
  Chi_K[ , pval_inform := pchisq( chisq_inform_std, df = 1L ) ]
  Chi_K[ , pval_uninform := pchisq( chisq_uninform_std, df = 1L ) ]
  Chi_K[ , `:=` ( df = NULL, chisq_inform = NULL, chisq_uninform = NULL  ) ]
  Chi_K[ , Fstat := chisq_inform_std / chisq_uninform_std ]
  Chi_K[ , pval_Fstat := pf( Fstat, df1 = 1, df2 = 1 ) ]
  Chi_K[]
  
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  
  if ( !smoothed ) {
    prop_H0_chisq = nrow( Chi_K[ pval_inform >= alpha ] ) / nrow( Chi_K )
    if ( prop_H0_chisq == 1 ) {
      cat( "---\n" )
      cat( "Null hypothesis is always accepted at a level of", alpha, "\n" )
    } else {
      reject = Chi_K[ pval_inform < alpha ]
      cat( "---\n" )
      cat( "Null hypothesis is rejected at a level of", alpha, 
           "for the following models\n" )
      cat( "---\n" )
      print( reject )
      cat( "---\n" )
      cat( "Overall, aggregated document stability is achieved for ", 
           round( prop_H0_chisq * 100, 2 ), "% of the models\n", sep = "" )
    }
  } else {
    smoothing = loess(data = Chi_K, chisq_inform_std ~ topic)
    hat = fitted.values( smoothing )
    smooth_dt = data.table( topic = Chi_K$topic, 
                            id_doc = Chi_K$id_doc, 
                            chisq_smooth = hat )
    smooth_dt[ , pval_inform := pchisq( chisq_smooth, df = 1L ) ]
    smooth_test = smooth_dt[ , .( chisq_smooth = mean( chisq_smooth ),
                                  pval_inform = mean(pval_inform) ), 
                             by = topic ]
    prop_H0_chisq = nrow( smooth_test[ pval_inform >= alpha ] ) / nrow( smooth_test )
    if ( prop_H0_chisq == 1 ) {
      cat( "---\n" )
      cat( "Null hypothesis is always accepted at a level of", alpha, "\n" )
    } else {
      reject = smooth_test[ pval_inform < alpha ]
      cat( "---\n" )
      cat( "Null hypothesis is rejected at a level of", alpha, 
           "for the following models\n" )
      cat( "---\n" )
      print( reject )
      cat( "---\n" )
      cat( "Overall, aggregated document stability is achieved for ", 
           round( prop_H0_chisq * 100, 2 ), "% of the models\n", sep = "" )
    }
  }
  
  if ( do_plot ) {
    cat( "Plotting...\n" )
    if ( !smoothed ) {
      p1 = ggplot( Chi_K ) +
        geom_hline( yintercept = qchisq( alpha, 1L ), linetype = 4 ) +
        geom_line( aes( x = topic, y = chisq_inform_std, color = as.factor(id_doc) ) ) +
        xlab( "Topics" ) + ylab( expression( bold( chi^2 ) ) ) +
        ggtitle( "Point-wise Aggregated Document Stability Plot" ) +
        theme_OpTop +
        theme( legend.position = "none" )
      p2 = ggplot( Chi_K ) +
        geom_hline( yintercept = qf( alpha, df1 = 1L, df2 = 1L ), linetype = 4 ) +
        geom_line( aes( x = topic, y = Fstat, color = as.factor(id_doc) ) ) +
        xlab( "Topics" ) + ylab( "F-stat" ) +
        ggtitle( "Point-wise F-Test Plot" ) +
        theme_OpTop +
        theme( legend.position = "none" )
      print( marrangeGrob( list( p1, p2 ), ncol = 1, nrow = 2, top = "" ) )
      
    } else {
      p1 = ggplot( Chi_K ) +
        geom_hline( yintercept = pchisq( alpha, 1L ), linetype = 4 ) +
        geom_smooth( aes( x = topic, y = chisq_inform_std ) ) +
        xlab( "Topics" ) + ylab( expression( bold( chi^2 ) ) ) +
        ggtitle( "Smoothed Aggregated Document Stability Plot" ) +
        theme_OpTop +
        theme( legend.position = "none" )
      p2 = ggplot( Chi_K ) +
        geom_hline( yintercept = qf( alpha, df1 = 1L, df2 = 1L ), linetype = 4 ) +
        geom_smooth( aes( x = topic, y = Fstat ) ) +
        xlab( "Topics" ) + ylab( "F-stat" ) +
        ggtitle( "Smoothed F-Test Plot" ) +
        theme_OpTop +
        theme( legend.position = "none" )
      print( marrangeGrob( list( p1, p2 ), ncol = 1, nrow = 2, top = "" ) )
      
    }
  }
  
  toc <- proc.time()
  runtime <- toc - tic
  cat( "---\n" )
  cat( "Function took:", runtime[ 3L ], "sec.\n" )
  cat( "---\n" )
  return( Chi_K )
  
}
