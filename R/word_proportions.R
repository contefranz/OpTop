if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "stopwords", "document", "id_doc",
                             "id_word", "word" ) )
}
#' Compute word proportions from a corpus object
#'
#' Fast routine that computes word proportions from a
#' \code{\link[quanteda]{corpus}} or \code{\link[quanteda]{dfm}} object.
#'
#' @param x Either a \code{\link[quanteda]{corpus}} or 
#' \code{\link[quanteda]{dfm}} object as defined in \code{quanteda}.
#' @param remove_document Remove the \code{document} identification inherited
#' from \code{x}. Default to \code{TRUE}.
#' @param remove_nonASCII A logical to remove non-ASCII characters from the 
#' \code{x}. Default to \code{TRUE}.
#' @param ... When \code{x} is a corpus, additional arguments passed 
#' to \code{\link[quanteda]{tokens}} and \code{\link[quanteda]{dfm}} to allow
#' for precision in tokens removal.
#' @details The function only applies preprocessing when \code{x} is a 
#' \code{\link[quanteda]{corpus}}. You can pass the usual parameters to \code{...}.
#' This includes the two parameters \code{remove_document} and \code{remove_nonASCII}.
#' If \code{x} is \code{\link[quanteda]{dfm}}, then the function directly
#' computes the word proportions as expected. 
#' @return A \code{data.table} with the following columns: 
#' \item{\code{id_doc}}{A sequential integer giving the identification of documents.}
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

word_proportions = function( x, remove_document = TRUE, 
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
    # word_prop <- dfm_weight( mydfm, "count" )
    word_count <- dfm_weight( mydfm, "count" )
    word_prop <- dfm_weight( mydfm, "prop" )
  } else if ( is.dfm( x ) ) {
    # word_prop <- dfm_weight( x, "count" )
    word_count <- dfm_weight( x, "count" )
    word_prop <- dfm_weight( x, "prop" )
  }
  
  # word_prop <- convert( word_prop, "data.frame" )
  # setDT( word_prop )
  word_count <- convert( word_count, "data.frame" )
  word_prop <- convert( word_prop, "data.frame" )
  setDT( word_count )
  setDT( word_prop )
  
  # word_prop = word_prop
  temp_count <- melt( word_count, id.vars = "document" )
  setorder( temp_count, document )
  setnames( temp_count, "variable", "word" )
  setnames( temp_count, "value", "word_count" )
  temp_count[ , id_doc := .GRP, by = document ]
  temp_count[ , id_word := 1L:.N, by = document ]
  # word_prop <- dfm_weight( mydfm, "prop" )
  # word_prop <- convert( word_prop, "data.frame" )
  # setDT( word_prop )
  # word_prop <- word_prop
  temp_prop <- melt( word_prop, id.vars = "document" )
  setorder( temp_prop, document )
  setnames( temp_prop, "variable", "word" )
  setnames( temp_prop, "value", "word_prop" )
  setkey( temp_count, document, word )
  setkey( temp_prop, document, word )
  out <- temp_count[ temp_prop, nomatch = 0L ]
  
  if ( remove_document ) {
    out[ , document := NULL ]
    setcolorder( out, c( "id_doc", "id_word", "word",
                         "word_count", "word_prop" ) )
  } else {
    setcolorder( out, c( "document", "id_doc", "id_word", "word",
                         "word_count", "word_prop" ) )
  }
  setkey( out, id_doc )
  return( out[] )
}
