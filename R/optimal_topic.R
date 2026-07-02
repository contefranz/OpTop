if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "id_word", "id_doc", "weighted_dfm", 
                             "word_prop", "document",
                             "word_sum", "check", "topics",
                             ".", "chisquare", "chisquare_mod",
                             "row_cut", "chi_sum", "word_prop_hat",
                             "word_prop_hat_cum", "pval", "docid",
                             "OpTop") )
}
#' Find the optimal number of topics from a pool of LDA models
#'
#' Identify the number of topics that best describes the corpus by applying a
#' fast chi-square–style test across a grid of `topicmodels::LDA` fits. The
#' routine evaluates each model and selects the optimal topic count using a
#' significance rule controlled by `alpha`, with a fallback to the global
#' minimum of the statistic.
#' 
#' @param lda_models A list of `topicmodels::LDA` objects (VEM), ordered by
#'   increasing number of topics. The grid should span the candidate values of `K`.
#' @param weighted_dfm A weighted `quanteda::dfm` containing word proportions
#'   for each document; it is recommended that document ids are available via
#'   `quanteda::docid()`.
#' @param q Numeric in `(0, 1]`. Cumulative mass used to define the truncated
#'   envelope (default `0.80`).
#' @param alpha Numeric in `[0, 1]`. Significance level for the selection rule
#'   (default `0.05`). See Details.
#' @param do_plot Logical; if `TRUE`, plot the statistic versus topics with
#'   vertical and horizontal guides at the selected optimum (default `TRUE`).
#'   
#' @details
#' The method builds a Pearson-type statistic under a multinomial view of word
#' counts and evaluates stability for each candidate `K`. The envelope used in
#' the comparison is truncated at cumulative mass `q` (e.g., `q = 0.80` keeps
#' the top 80% of mass).
#'
#' **Selection rule.** The optimal `K` is the first model whose p-value is at
#' most `alpha`. If no model reaches `alpha`, the model with the minimum
#' statistic is selected. Setting `alpha = 0` forces the global-minimum rule.
#'
#' **Input alignment.** `weighted_dfm` must be a `quanteda::dfm` of word
#' proportions (row-wise). Document identifiers are taken from
#' `quanteda::docid(weighted_dfm)` and matched to the LDA fits; documents not
#' present in the first model’s `@documents` slot are dropped (with a message)
#' to ensure alignment. The remaining documents are paired with the rows of
#' `@gamma` by identifier, so the row order of `weighted_dfm` does not need to
#' follow the order in which the models saw the documents.
#'
#' **Performance note.** The core computation is delegated to C++ compiled code
#' to handle high-dimensional vocabularies efficiently.
#' 
#' @return A `data.table` with columns:
#' - `topic`: integer number of topics (`K`).
#' - `OpTop`: standardized chi-square–style statistic for each `K`.
#' - `pval`: p-value associated with `OpTop`.
#' 
#' @examples
#' \dontrun{
#' # Compute word proportions from a corpus objects
#' test1 = optimal_topic( lda_models = lda_list,
#'                         weighted_dfm = weighted_dfm,
#'                         q = 0.80,
#'                         alpha = 0.05 )
#' }
#' 
#' @seealso [topicmodels::LDA()]
#' 
#' @import data.table ggplot2
#' @importFrom quanteda is.dfm docid
#' @export

optimal_topic = function( lda_models, weighted_dfm, q = 0.80, alpha = 0.05, do_plot = TRUE ) {
  
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
  
  tic = proc.time()
  docs = as.character( docid( weighted_dfm ) )

  # drop documents the models never saw; all models are assumed to share the
  # document set of the first one (a document the LDA drops is dropped for
  # every k), so the check runs against lda_models[[ 1L ]] only
  doc_check = docs %in% lda_models[[ 1L ]]@documents
  if ( !all( doc_check ) ) {
    id_toremove = which( !doc_check )
    if ( length( id_toremove ) < length( doc_check ) ) {
      cat("Removing unmatched documents\n" )
      weighted_dfm = weighted_dfm[ -id_toremove, ]
      docs = docs[ -id_toremove ]
    } else {
      stop("Document matching went really wrong. Check docs in both weighted_dfm and in LDA@documents")
    }
  }

  # map each dfm row to the corresponding row of @gamma (0-based for C++);
  # membership was checked above, so no NA can survive the match
  doc_map = match( docs, lda_models[[ 1L ]]@documents ) - 1L

  # cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Beginning computations...\n" )
  Chi_K = optimal_topic_core( lda_models, weighted_dfm, q, doc_map )
  # cat( "# # # # # # # # # # # # # # # # # # # #\n" )
  cat( "Computations done!\n" )
  cat( "---\n" )
  
  Chi_K = as.data.table(Chi_K)
  setnames( Chi_K, old = names( Chi_K ), c( "topic", "OpTop", "pval" ) )  
  
  global_min = Chi_K[ , .SD[ which.min( OpTop ) ] ]
  alpha_min = Chi_K[ pval <= alpha ][ 1L ]
  if ( alpha == 0 || all( is.na( alpha_min ) ) ) {
    cat( "Optimal model found by global minimum\n" )
    cat( "Optimal model has", global_min$topic, "topics\n" )
    best_topic = global_min
  } else {
    if ( global_min$topic > alpha_min$topic ) {
      cat( "Optimal model found by significance level of", alpha, "\n" )
      cat( "Optimal model has", alpha_min$topic, "topics\n" )
      best_topic = alpha_min
    } else {
      cat( "Optimal model found by global minimum\n" )
      cat( "Optimal model has", global_min$topic, "topics\n" )
      best_topic = global_min
    }
  }
  
  if ( do_plot ) {
    cat( "Plotting...\n" )
    x_min = best_topic$topic
    y_min = best_topic$OpTop
    p1 = ggplot( Chi_K ) +
      geom_line( aes( x = topic, y = OpTop ), linewidth = 0.8, color = "royalblue" ) +
      geom_hline( yintercept = y_min, color = "black", linetype = 2L ) +
      geom_vline( xintercept = x_min, color = "black", linetype = 2L ) +
      annotate( "point", x = x_min, y = y_min, color = "red", shape = 4L, size = 4L ) +
      xlab( "Topics" ) + ylab( expression(OpTop[J]^{"K"}) ) +
      ggtitle( "Optimal Topic Plot" ) +
      theme_OpTop
    print( p1 )
  }
  
  toc = proc.time()
  runtime = toc - tic
  cat( "---\n" )
  cat( "Function took:", runtime[ 3L ], "sec.\n" )
  cat( "---\n" )
  return( Chi_K )
}
