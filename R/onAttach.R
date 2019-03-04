#' @importFrom utils packageVersion
.onAttach = function( libname, pkgname ) {
  # Runs when attached to search() path such as by library() or require()
  if ( interactive() ) {
    packageStartupMessage( 'Welcome to OpTop version ',
                           as.character( packageVersion( "OpTop" ) ) )
    # packageStartupMessage( 'For help ?OpTop or vignette( "OpTop" )' )
    packageStartupMessage( 'If you find bugs, please report them at https://github.com/contefranz/OpTop/issues' )
  }
}
