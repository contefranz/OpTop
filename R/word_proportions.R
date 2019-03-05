if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables( c( "stopwords", "document", "id_doc",
                             "id_word", "word" ) )
}
#' Compute word proportions from a corpus object
#'
#' Fast routine that computes word proportions from a
#' \code{\link[quanteda]{corpus}}.
#'
#' @param corpus A \code{\link[quanteda]{corpus}} object as defined in
#' \code{quanteda}.
#' @param remove_document Remove the \code{document} identification inherited
#' from \code{\link[quanteda]{corpus}}. Default to \code{TRUE}.
#' @param language A character specifying the language of stopwords according
#' to ISO 639-1 code. It mimics the behavior of \code{\link[stopwords]{stopwords}}.
#' Default to \code{"en"} for English.
#' @param source A character specifying the stopwords source. It can be either
#' \code{"snowball"}, \code{"snowball-iso"}, \code{"misc"}, or \code{"smart"}.
#' Default to \code{"snowball"} for consistency with
#' \code{\link[stopwords]{stopwords}}.
#' @return A \code{data.table} with the following columns:
#' \item{\code{id_doc}}{A sequential integer giving the identification of documents.}
#' \item{\code{id_word}}{A sequential integer giving the identification of words.}
#' \item{\code{word}}{A character identifying the word.}
#' \item{\code{word_count}}{An integer giving the word count.}
#' \item{\code{word_prop}}{A numeric giving the word proportion.}
#' @examples
#' \dontrun{
#' # Compute word proportions from a corpus objects
#' word_proportions = word_proportions( corpus = data_corpus_inaugural,
#'                                      remove_document = TRUE,
#'                                      language = "en",
#'                                      source = "snowball" )
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
#' @importFrom quanteda dfm
#' @importFrom quanteda stopwords
#' @importFrom quanteda dfm_weight
#' @importFrom quanteda convert
#' @export

word_proportions = function( corpus, remove_document = FALSE,
                             language = "en", source = "snowball" ) {

  if ( !is.corpus( corpus ) ) {
    stop( "corpus must be a corpus object as defined by quanteda")
  }
  if ( !is.logical( remove_document ) ) {
    stop( "remove_document must be either TRUE or FALSE" )
  }
  if ( !is.character( language ) ) {
    stop( "language must a character" )
  }
  if ( !is.character( source ) ) {
    stop( "source must be a character" )
  }

  mydfm = dfm( x = corpus,
               remove = stopwords( language = language, source = source ),
               remove_punct = TRUE,
               remove_numbers = TRUE,
               remove_symbols = TRUE )

  word_prop = dfm_weight( mydfm, "count" )
  word_prop = convert( word_prop, "data.frame" )
  setDT( word_prop )

  word_prop = word_prop
  temp1 = melt( word_prop, id.vars = "document" )
  setorder( temp1, document )
  setnames( temp1, "variable", "word" )
  setnames( temp1, "value", "word_count" )
  temp1[ , id_doc := .GRP, by = document ]
  temp1[ , id_word := 1L:.N, by = document ]
  word_prop = dfm_weight( mydfm, "prop" )
  word_prop = convert( word_prop, "data.frame" )
  setDT( word_prop )
  word_prop = word_prop
  temp2 = melt( word_prop, id.vars = "document" )
  setorder( temp2, document )
  setnames( temp2, "variable", "word" )
  setnames( temp2, "value", "word_prop" )
  setkey( temp1, document, word )
  setkey( temp2, document, word )
  out = temp1[ temp2, nomatch = 0L ]

  if ( remove_document ) {
    out[ , document := NULL ]
    setcolorder( out, c( "id_doc", "id_word", "word",
                         "word_count", "word_prop" ) )
  } else {
    setcolorder( out, c( "document", "id_doc", "id_word", "word",
                         "word_count", "word_prop" ) )
  }
  return( out[] )
}
