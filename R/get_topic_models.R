#' Get the list of topic models from a specified environment
#'
#' Easily get the list of topic models from a specific class. This is a handy
#' feature to facilitate the identification of only the objects the user wants
#' to process.
#' @inheritParams base::ls
#' @param object_type A character specifying the class of objects to extract.
#' Default to `LDA_VEM` as given by the virtual class
#' [topicmodels::TopicModel-class()].
#' @param envir The environment where to search for topic models. Default to
#' [base::globalenv()].
#' @return A named list containing the objects of specified class.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}.
#' @export

get_topic_models <- function( object_type, pattern, envir = globalenv() ) {

  if( !missing( object_type ) && !is.character( object_type ) ) {
    stop( "object_type must be a character" )
  }
  if ( !missing( pattern ) && !is.character( pattern ) ) {
    stop( "pattern must be a regular expression given as character" )
  }
  if ( !is.environment( envir ) ) {
    stop( "envir must be an environment where to search for topic models" )
  }
  if( missing( pattern ) ) {
    pattern = ""
  }

  obj = objects( pattern = pattern, envir = envir )
  # if ( missing( object_type ) ) {
  #   # pos = sapply( ls( pattern = pattern, envir = envir ),
  #   #               function( x ) is.LDA_VEM( get( x ) ) )
  # } else {
  #   pos = sapply( ls( pattern = pattern, envir = envir ),
  #                 function( x ) class( get( x ) ) == object_type )
  # }
  # 
  # models_names = ls( envir = envir )[ pos ]
  models_list = mget( x = obj, envir = envir )
  # if ( missing( object_type ) ) {
  #   all( sapply( models_list, is.LDA_VEM ) )
  # }
  
  return( models_list )
}


