#' Cumulative sum scaling factors.
#' 
#' Calculates each column's quantile and calculates the sum up to and including
#' that quantile.
#' 
#' 
#' @param obj A MRexperiment object.
#' @param p The pth quantile.
#' @return Returns a matrix normalized by scaling counts up to and including
#' the pth quantile.
#' @seealso \code{\link{fitZig}} \code{\link{cumNorm}}
#' @examples
#' 
#' data(mouseData)
#' head(cumNormMat(mouseData))
#' 
cumNormMat <-
function(obj,p= cumNormStat(obj)){
####################################################################################
#   Calculates each column's quantile
#    and calculated the sum up to and
#    including that quantile.
####################################################################################
    x=MRcounts(obj)
    xx=x
	xx[x==0] <- NA
	
	qs=matrixStats::colQuantiles(xx,p=p,na.rm=TRUE)
	
	newMat<-sapply(1:ncol(xx), function(i) {
				   xx=(x[,i]-.Machine$double.eps)
				   sum(xx[xx<=qs[i]])
				   })
	x<-sweep(x,2,newMat/1000,"/")
	return(x)
}
