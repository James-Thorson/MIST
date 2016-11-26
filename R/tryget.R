
#' Try to get an object from memory
#'
#' \code{tryget} checks for an object and returns it if its available

tryget = function( char ){
  if( exists(char) ){
    Return = get(char)
  }else{
    Return = NULL
  }
  return( Return )
}
