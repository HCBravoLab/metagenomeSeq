#' Cumulative sum scaling factors.
#'
#' Calculates each column's quantile and calculates the sum up to and including that quantile.
#'
#' @param jobj An eSet object.
#' @param p The pth quantile.
#' @return jobj An updated eSet object with normalized counts
#'
#' @name cumNormMat
#' @seealso \code{\link{fitZig}} \code{\link{cumNorm}}
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
