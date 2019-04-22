#' Cumulative sum scaling normalization
#' 
#' Calculates each column's quantile and calculates the sum up to and including
#' that quantile.
#' 
#' @param obj An MRexperiment object.
#' @param p The pth quantile.
#' @return Object with the normalization factors stored as 
#' a vector of the sum up to and including a sample's pth quantile.
#' @seealso \code{\link{fitZig}} \code{\link{cumNormStat}}
#' @examples
#' 
#' data(mouseData)
#' mouseData <- cumNorm(mouseData)
#' head(normFactors(mouseData))
#' 
cumNorm <- function(obj,p=cumNormStatFast(obj)){
	if(class(obj)=="MRexperiment"){
		x = MRcounts(obj,norm=FALSE,log=FALSE)
	} else {
		stop("Object needs to be a MRexperiment object")
	}
	normFactors = calcNormFactors(obj=x,p=p)
	pData(obj@expSummary$expSummary)$normFactors = normFactors
	validObject(obj)
	return(obj)
}

#' Cumulative sum scaling (css) normalization factors
#' 
#' Return a vector of the the sum up to and including a quantile.
#' 
#' @param obj An MRexperiment object or matrix.
#' @param p The pth quantile.
#' @return Vector of the sum up to and including a sample's pth quantile.
#' @seealso \code{\link{fitZig}} \code{\link{cumNormStatFast}} \code{\link{cumNorm}}
#' @examples
#' 
#' data(mouseData)
#' head(calcNormFactors(mouseData))
#' 
calcNormFactors <- function(obj,p=cumNormStatFast(obj)){
	x = returnAppropriateObj(obj,norm=FALSE,log=FALSE)
	xx = x
	xx[x == 0] <- NA
	qs = colQuantiles(xx, probs = p, na.rm = TRUE)
	normFactors <- sapply(1:ncol(xx), function(i) {
		xx = (x[, i] - .Machine$double.eps)
		sum(xx[xx <= qs[i]])
	})
	names(normFactors)<-colnames(x)
	as.data.frame(normFactors)
}
