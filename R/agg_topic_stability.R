if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "chisq_std", "df", ".") )
}
#' Compute Aggregate Topic Stability
#'
#' `r lifecycle::badge("deprecated")`
#'
#' As of OpTop 0.9.8, `agg_topic_stability()` is deprecated and scheduled for
#' removal: the package is converging on the discrepancy-index API (see
#' [optop_index_table()]).
#'
#' Compute aggregate topic stability and perform a chi-square test to evaluate
#' how topic structure changes for models with more topics than the selected
#' optimal model. The function produces a standardized chi-square statistic
#' per document and model, and optional plots showing how stability varies with
#' the number of topics.
#' 
#' @param lda_models A list of fitted `topicmodels::LDA` models (VEM), typically
#'   ordered by increasing number of topics and including `optimal_model`.
#' @param optimal_model Integer giving the number of topics for the optimal model;
#'   must match one entry in `lda_models`.
#' @param q Numeric in `(0, 1]`; cumulative mass used to define the best-pair
#'   envelope (default `0.80`).
#' @param alpha Numeric significance level used when assessing stability
#'   (default `0.05`).
#' @param smoothed Logical; if `TRUE`, apply LOESS smoothing across documents
#'   before testing (default `TRUE`). If `FALSE`, use pointwise document-level
#'   values.
#' @param do_plot Logical; if `TRUE`, print the chi-square plot as a function of
#'   the number of topics (default `TRUE`).
#'
#' @details
#' For each model with `k > optimal_model`, the routine compares the best-matching
#' topics to the reference structure at `optimal_model` and builds, for each
#' document, a cumulative “best-pair” envelope up to mass `q`. The standardized
#' chi-square statistic is then computed on that envelope, and a p-value is
#' obtained from the chi-square distribution with one degree of freedom.
#'
#' When `smoothed = TRUE`, a LOESS smoother is fit to the document-level
#' chi-square values across topics; inference is then based on the smoothed
#' series. When `smoothed = FALSE`, inference is based on the pointwise
#' document-level values. If `do_plot = TRUE`, the function prints a plot of
#' the standardized chi-square against the number of topics (with either
#' per-document trajectories or a smoothed curve, depending on `smoothed`).
#'
#' Inputs are expected to be fitted [topicmodels::LDA] objects (VEM). The
#' `optimal_model` must match one of the topic counts present in `lda_models`.
#' If the optimal model is already the last element (largest k) in `lda_models`,
#' there is nothing to evaluate above it and the function returns `NULL`
#' with a message.
#' 
#' @return A `data.table` with one row per document–model pair and the columns:
#' - `topic`: integer number of topics.
#' - `id_doc`: integer document id (row index in the aligned data).
#' - `chisq_std`: standardized chi-square statistic.
#' - `pval`: p-value of the chi-square test.
#' 
#' @examples
#'\dontrun{
#' test4 <- agg_topic_stability( lda_models = lda_list,
#'                               optimal_model = test1)
#' }
#' 
#' @seealso [optimal_topic()] [topicmodels::LDA()]
#' 
#' @import data.table
#' @export

agg_topic_stability <- function( lda_models, optimal_model,
                                 q = 0.80, alpha = 0.05,
                                 smoothed = TRUE,
                                 do_plot = TRUE ) {

  lifecycle::deprecate_warn(
    when = "0.9.8", what = "agg_topic_stability()",
    details = paste( "OpTop is converging on the discrepancy-index API",
                     "(see optop_index_table()); this function will be",
                     "removed in a future release." )
  )
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
  if ( !is.numeric( q ) ) {
    stop( "q must be a numeric" )
  }
  if ( !is.numeric( alpha ) ) {
    stop( "alpha must be a numeric" )
  }
  if ( !is.logical( smoothed ) ) {
    stop( "smoothed must be either TRUE or FALSE" )
  }
  if ( optimal_model == lda_models[[ length( lda_models ) ]]@k ) {
    message("Optimal model is already the last one in lda_models. There is nothing to compute above that.")
    return( NULL )
  }
  
  ###############################################################
  # THIS IS OLD CODE TO BE REMOVED ONCE THE CONVERSION IS STABLE
  ###############################################################
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
  ###############################################################
  
  ###############################################################
  # THIS IS OLD CODE TO BE REMOVED ONCE THE CONVERSION IS STABLE
  ###############################################################
  # cat( "best model has", .optimal_model, "topics\n" )
  # tic <- proc.time()
  # 
  # k_end <- max( sapply( lda_models, function( x ) x@k ) )
  # best_pos <- which( sapply( lda_models, function( x ) x@k ) == .optimal_model )
  # 
  # if ( length( best_pos ) == 0 ) {
  #   stop( paste( "There is no optimal model in lda_models.",
  #                "This could be either due to a wrong specification of",
  #                "argument optimal_model or",
  #                "if optimal_model is a data.table, the optimal model cannot be found",
  #                "in the list lda_models." ) )
  # }
  ###############################################################
  tic <- proc.time()
  best_pos <- which( sapply( lda_models, function( x ) x@k ) == optimal_model )
  if ( length( best_pos ) == 0 ) {
    stop("optimal_model does not correspond to any topic number in lda_models")
  }
  
  ##########################
  # C++ BEGINS HERE ! ! !
  ##########################
  
  if ( !is.LDA_VEM( optimal_model ) ) {
    # extracting information from best model
    dtw_best <- lda_models[[ best_pos ]]@gamma
    tww_best <- t( exp( lda_models[[ best_pos ]]@beta ) )
  }
  n_doc <- nrow( dtw_best )
  loop_sequence <- (best_pos + 1L):length( lda_models )  
  min_loop <- min( loop_sequence )
  k = 1L
  Chi_K <- data.table()
  cat( "---\n" )
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  matches = topic_match(lda_models = lda_models, optimal_model = optimal_model)

  for ( i_mod in loop_sequence ) {
    # i_pos <- i_mod - .optimal_model + 1L
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
    # i <- current_k - .optimal_model
    MostSim <- intersect( 1L:current_k, matches$BestMatch[ i, ] )
    MostSimDTW <- dtw[ , MostSim ]
    MostSimTWW <- tww[ , MostSim ]
    
    out_in <- data.table()
    cat( "--> Processing documents\n" )
    for ( j_doc in 1L:n_doc ) {
      DTW_best <- matrix( dtw_best[ j_doc, ], 
                          nrow = n_tww,
                          ncol = optimal_model )
      DTW_most <- matrix( MostSimDTW[ k, ], 
                          nrow = n_tww,
                          ncol = ncol( MostSimDTW ) )
      
      tww_dtw_best <- DTW_best * tww_best
      ones_tww_dtw_best <- matrix( 1., nrow = ncol( tww_dtw_best ), ncol = 1L )
      totprob_best <- sum( tww_dtw_best %*% ones_tww_dtw_best, na.rm = TRUE )
      tww_dtw_best2 <- tww_dtw_best %*% rep( 1L, ncol( tww_dtw_best ) ) / totprob_best
      
      tww_dtw <- DTW_most * MostSimTWW
      totprob <- sum( tww_dtw %*% rep( 1L, ncol( tww_dtw ) ) )
      tww_dtw2 <- tww_dtw %*% rep( 1L, ncol( tww_dtw ) ) / totprob
      
      X <- cbind( tww_dtw_best2, tww_dtw2 )
      BestPair <- apply( X, 2L, function( x ) sort( x, decreasing = TRUE ) )
      BestPair <- cbind( BestPair, cumsum( BestPair[ , 1L ] ) )
      n_BP <- nrow( BestPair )
      p_BP <- ncol( BestPair )
      # stop when you reach q
      AggBestPair <- BestPair[ BestPair[ , 3L ] <= q, ]
      icut <- nrow( AggBestPair )
      
      # icut <- base::which.min( abs( BestPair[ , 1L ] - q ) )
      # if ( icut > 250 ) {
      lowest_estimates <- apply( BestPair[ (icut + 1L):n_BP, ], 2L, sum )
      # sum_overbest <- apply( BestPair[ (icut + 1L):nrow( BestPair ), ], 2L, sum )
      AggBestPair <- rbind( BestPair[ 1L:icut, ], lowest_estimates )
      
      numerator <- ( AggBestPair[ , 1L ] - AggBestPair[ , 2L ] )^2L
      denominator <- AggBestPair[ , 1L ]
      chisq <- icut * sum( numerator / denominator, na.rm = TRUE )
      
      out_doc <- data.table( topic = current_k,
                             id_doc = j_doc,
                             df = icut,
                             chisq = chisq )
      out_in <- rbindlist( list( out_in, out_doc ) )
      # }
    }
    Chi_K <- rbindlist( list( Chi_K, out_in ) )
  }
  Chi_K[ , chisq_std := chisq / df ]
  Chi_K[ , pval := pchisq( chisq_std, df = 1L ) ]
  Chi_K[ , `:=` ( df = NULL, chisq = NULL ) ]
  Chi_K[]
  
  cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  
  if ( !smoothed ) {
    prop_H0 = nrow( Chi_K[ pval >= alpha ] ) / nrow( Chi_K )
    if ( prop_H0 == 1 ) {
      cat( "---\n" )
      cat( "Null hypothesis is always accepted at a level of", alpha, "\n" )
    } else {
      reject = Chi_K[ pval < alpha ]
      cat( "---\n" )
      cat( "Null hypothesis is rejected at a level of", alpha, 
           "for the following models\n" )
      cat( "---\n" )
      print( reject )
      cat( "---\n" )
      cat( "Overall, aggregated topic stability is achieved for ", 
           round( prop_H0 * 100, 2 ), "% of the models\n", sep = "" )
    }
  } else {
    smoothing = loess(data = Chi_K, chisq_std ~ topic)
    hat = fitted.values( smoothing )
    smooth_dt = data.table( topic = Chi_K$topic, 
                            id_doc = Chi_K$id_doc, 
                            chisq_smooth = hat )
    smooth_dt[ , pval := pchisq( chisq_smooth, df = 1) ]
    smooth_test = smooth_dt[ , .( chisq_smooth = mean( chisq_smooth ),
                                  pval = mean(pval) ), 
                             by = topic ]
    prop_H0 = nrow( smooth_test[ pval >= alpha ] ) / nrow( smooth_test )
    if ( prop_H0 == 1 ) {
      cat( "---\n" )
      cat( "Null hypothesis is always accepted at a level of", alpha, "\n" )
    } else {
      reject = smooth_test[ pval < alpha ]
      cat( "---\n" )
      cat( "Null hypothesis is rejected at a level of", alpha, 
           "for the following models\n" )
      cat( "---\n" )
      print( reject )
      cat( "---\n" )
      cat( "Overall, aggregated topic stability is achieved for ", 
           round( prop_H0 * 100, 2 ), "% of the models\n", sep = "" )
    }
  }
  
  if ( do_plot ) {
    cat( "Plotting...\n" )
    if ( !smoothed ) {
    p1 = ggplot( Chi_K ) +
      geom_hline( yintercept = qchisq( alpha, 1L ), linetype = 4 ) +
      geom_line( aes( x = topic, y = chisq_std, color = as.factor(id_doc) ) ) +
      xlab( "Topics" ) + ylab( expression( bold( chi^2 ) ) ) +
      ggtitle( "Point-wise Aggregated Topic Stability Plot" ) +
      theme_OpTop +
      theme( legend.position = "none" )
    print( p1 )
    } else {
      p1 = ggplot( Chi_K ) +
        geom_hline( yintercept = pchisq( alpha, 1L ), linetype = 4 ) +
        geom_smooth( aes( x = topic, y = chisq_std ) ) +
        xlab( "Topics" ) + ylab( expression( bold( chi^2 ) ) ) +
        ggtitle( "Smoothed Aggregated Topic Stability Plot" ) +
        theme_OpTop +
        theme( legend.position = "none" )
      print( p1 )
    }
  }
  
  toc <- proc.time()
  runtime <- toc - tic
  cat( "---\n" )
  cat( "Function took:", runtime[ 3L ], "sec.\n" )
  cat( "---\n" )
  return( Chi_K )
}
