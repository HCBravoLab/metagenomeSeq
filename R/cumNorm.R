#' Cumulative sum scaling factors.
#' 
#' Calculates each column's quantile and calculates the sum up to and including
#' that quantile.
#' 
#' 
#' @param obj An MRexperiment object.
#' @param p The pth quantile.
#' @return Vector of the sum up to and including a sample's pth quantile
#' @seealso \code{\link{fitZig}} \code{\link{cumNormStat}}
#' @examples
#' 
#' data(mouseData)
#' cumNorm(mouseData)
#' head(normFactors(mouseData))
#' 
cumNorm <-
function(obj,p=cumNormStat(obj)){
	x = MRcounts(obj,norm=FALSE,log=FALSE)
	xx=x
	xx[x==0] <- NA
		
	qs=matrixStats::colQuantiles(xx,p=p,na.rm=TRUE)
		
	normFactors<-sapply(1:ncol(xx), function(i) {
		   xx=(x[,i]-.Machine$double.eps)
		   sum(xx[xx<=qs[i]])
		   })
	names(normFactors)<- colnames(x)
        pData(obj@expSummary$expSummary)$normFactors = as.data.frame(normFactors)
        validObject(obj)
	return(obj)
}
