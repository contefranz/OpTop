# handy function to check LDA_VEM class
is.LDA_VEM <- function( x ) {
  "LDA_VEM" %in% class(x)
}
#' @keywords internal

# this declare my personal ggplot theme
font_size = 10
theme_OpTop <- theme( title = element_text( face = "bold", size = 8 ),
                      axis.title.x = element_text( face = "bold", size = font_size ),
                      axis.title.y = element_text( face = "bold", size = font_size ),
                      axis.text.x = element_text( size = font_size ),
                      axis.text.y = element_text( size = font_size ),
                      legend.text = element_text( size = 8 ),
                      legend.position = "bottom" )
