fixdiag = function( Ncol, Nrow ){
  tmp = matrix( 1:(Ncol*Nrow), nrow=Nrow, ncol=Ncol)
  diag(tmp) = NA
  tmp = tmp[lower.tri(tmp,diag=TRUE)]
  return( tmp )
}
