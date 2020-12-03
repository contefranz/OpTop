if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "stopwords", "doc_id", "id_doc",
                             "id_word", "word" ) )
}
#' Compute word proportions from a corpus or a dfm object
#'
#' Fast routine that computes word proportions from a
#' \code{\link[quanteda]{corpus}} or \code{\link[quanteda]{dfm}} object.
#'
#' @param x Either a \code{\link[quanteda]{corpus}} or 
#' \code{\link[quanteda]{dfm}} object as defined in \code{quanteda}.
#' @param remove_document Remove the \code{document} identification inherited
#' from \code{x}. Default to \code{FALSE} (See 'Details').
#' @param remove_nonASCII A logical to remove non-ASCII characters from the 
#' \code{x}. Default to \code{TRUE}.
#' @param ... When \code{x} is a corpus, additional arguments passed 
#' to \code{\link[quanteda]{tokens}} and \code{\link[quanteda]{dfm}} to allow
#' for precision in tokens removal.
#' @details We recommend to keep \code{remove_document = FALSE} since
#' this triggers further controls before estimating the optimal model specification.
#' We refer the user to the inline help of \code{\link[OpTop]{optimal_topic}}.
#' 
#' \code{word_proportions} only applies preprocessing when \code{x} is a 
#' \code{\link[quanteda]{corpus}}. You can pass the usual parameters to \code{...}.
#' This includes the two parameters \code{remove_document} and \code{remove_nonASCII}.
#' If \code{x} is \code{\link[quanteda]{dfm}}, then the function directly
#' computes the word proportions as expected. 
#' @return A \code{data.table} with the following columns:
#' \item{\code{document}}{A character giving the original document identification
#' as inherited from the corpus or the dfm. This column is only retained when 
#' \code{remove_document = FALSE}.} 
#' \item{\code{id_doc}}{A sequential integer giving the document identification.}
#' \item{\code{id_word}}{A sequential integer giving the identification of words.}
#' \item{\code{word}}{A character identifying the word.}
#' \item{\code{word_count}}{An integer giving the word count.}
#' \item{\code{word_prop}}{A numeric giving the word proportion.}
#' @examples
#' \dontrun{
#' # Compute word proportions from a corpus object
#' word_proportions <- word_proportions( x = data_corpus_inaugural,
#'                                       remove = stopwords(), tolower = TRUE )
#' # Compute word proportions from a dfm object
#' word_proportions <- word_proportions( x = data_dfm_lbgexample )
#' }
#' @seealso \code{\link[quanteda]{corpus}} \code{\link[data.table]{data.table}}
#' \code{\link[stopwords]{stopwords}}
#' @references Stopwords ISO: \url{https://github.com/stopwords-iso/stopwords-iso}.\cr
#'
#' Full list of stopwords:
#' \url{https://github.com/stopwords-iso/stopwords-iso/blob/master/CREDITS.md}.
#' @author Francesco Grossetti \email{francesco.grossetti@@unibocconi.it}.
#' @import data.table
#' @importFrom quanteda is.corpus 
#' @importFrom quanteda dfm is.dfm
#' @importFrom quanteda stopwords
#' @importFrom quanteda dfm_weight
#' @importFrom quanteda convert
#' @export

word_proportions = function( x, remove_document = FALSE, 
                             remove_nonASCII = TRUE, ... ) {

  if ( !is.corpus( x ) && !is.dfm( x ) ) {
    stop( "x must be a corpus or a dfm object as defined by quanteda")
  }
  if ( !is.logical( remove_document ) ) {
    stop( "remove_document must be either TRUE or FALSE" )
  }
  if ( !is.logical( remove_nonASCII ) ) {
    stop( "remove_nonASCII must be either TRUE or FALSE" )
  }
  
  if ( is.corpus( x ) ) {
    mydfm <- dfm( x, ... )
    if ( remove_nonASCII ) {
      mydfm <- dfm( mydfm, remove = "[^ -~]", valuetype = "regex" )
    }
    word_count <- dfm_weight( mydfm, "count" )
    word_prop <- dfm_weight( mydfm, "prop" )
  } else if ( is.dfm( x ) ) {
    word_count <- dfm_weight( x, "count" )
    word_prop <- dfm_weight( x, "prop" )
  }
  
  word_count <- convert( word_count, "data.frame" )
  word_prop <- convert( word_prop, "data.frame" )
  setDT( word_count )
  setDT( word_prop )
  
  temp_count <- melt( word_count, id.vars = "doc_id", 
                      variable.name = "word", value.name = "word_count" )
  setorder( temp_count, doc_id )
  temp_count[ , id_doc := .GRP, by = doc_id ]
  temp_count[ , id_word := 1L:.N, by = doc_id ]
  temp_prop <- melt( word_prop, id.vars = "doc_id",
                     variable.name = "word", value.name = "word_prop")
  setorder( temp_prop, doc_id )
  setkey( temp_count, doc_id, word )
  setkey( temp_prop, doc_id, word )
  out <- temp_count[ temp_prop, nomatch = 0L ]
  
  if ( remove_document ) {
    out[ , doc_id := NULL ]
    setcolorder( out, c( "id_doc", "id_word", "word",
                         "word_count", "word_prop" ) )
  } else {
    setcolorder( out, c( "id_doc", "doc_id", "id_word", "word",
                         "word_count", "word_prop" ) )
  }
  setkey( out, NULL )
  return( out[] )
}
