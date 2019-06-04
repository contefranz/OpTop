# handy function to check LDA_VEM class
is.LDA_VEM <- function( x ) {
  "LDA_VEM" %in% class(x)
}

# this normalize the tww matrix corresponding to the optimal topic
# specification given by optimal_topic()
norm_tww <- function( tww ) {
  
  # Normalizing by scaling by vector norms
  norms <- apply( tww, 2L, function( x ) sqrt( sum( abs( x )^2L ) ) )
  tww_norm <- scale( tww, center = FALSE, scale = norms )
  return( tww_norm )
  
}

# this declare my personal ggplot theme
font_size = 10
theme_OpTop <- theme( title = element_text( face = "bold", size = 8 ),
                      axis.title.x = element_text( face = "bold", size = font_size ),
                      axis.title.y = element_text( face = "bold", size = font_size ),
                      axis.text.x = element_text( size = font_size ),
                      axis.text.y = element_text( size = font_size ),
                      legend.text = element_text( size = 8 ),
                      legend.position = "bottom" )
