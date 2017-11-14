
#' Plot interaction matrix
#'
#' \code{plot_interactions} plots the estimated interaction matrix
#'
#' @export
plot_interactions = function( Report, TmbData, category_names=1:TmbData$n_p, plotdir=paste0(getwd(),"/"), res=200, ... ){

  # Preparatory calculations
  Eigen = eigen(Report$B_pp)
  isComplexTF = ifelse( any(abs(Im(Eigen$values))>0.001), TRUE, FALSE)
  Ncol = ceiling(sqrt(TmbData$n_p))
  Nrow = ceiling(TmbData$n_p/Ncol)

  # Plot eigen-decomposition
  if( !is.null(plotdir)) png( file=paste0(plotdir,"eigen_B_pp.png"), width=3*Ncol, height=3*Nrow, units="in", res=res)
    par( mfrow=c(Nrow,Ncol), ... )
    for(p in 1:TmbData$n_p){
      if( isComplexTF==TRUE ){
        plot( x=1:TmbData$n_p-0.05, y=Re(Eigen$vectors[,p])*sign(mean(Re(Eigen$vectors[,p]))), type="h", lwd=3, ylim=c(-1,1.2), xaxt="n", xlab="", ylab="", xlim=c(1,TmbData$n_p)+c(-0.5,0.5) )
        lines( x=1:TmbData$n_p+0.05, y=Im(Eigen$vectors[,p])*sign(mean(Re(Eigen$vectors[,p]))), type="h", lwd=3, col="red")
        legend( "top", bty="n", legend=paste0(ThorsonUtilities::FormatC(Re(Eigen$values[p]),3)," ",ifelse(sign(Im(Eigen$values[p]))<0,"-","+")," ",FormatC(abs(Im(Eigen$values[p])),3),"i"), cex=1.2 )
      }else{
        plot( x=1:TmbData$n_p, y=Eigen$vectors[,p]*sign(mean(Eigen$vectors[,p])), type="h", lwd=3, ylim=c(-1,1.2), xaxt="n", xlab="", ylab="", xlim=c(1,TmbData$n_p)+c(-0.5,0.5) )
        legend( "top", bty="n", legend=FormatC(Eigen$values[p],3), cex=1.2 )
      }
      abline( h=0, lty="dotted")
      if( p > (TmbData$n_p-Ncol) ){
        axis(side=1, at=1:TmbData$n_p, line=0, labels=rep("",TmbData$n_p) )
        mtext(side=1, line=1.5, at=1:TmbData$n_p, text=category_names )
      }
    }
  if( !is.null(plotdir)) dev.off()
}
