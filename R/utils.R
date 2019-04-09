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
