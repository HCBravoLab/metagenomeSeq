#' Cumulative sum scaling percentile selection
#' 
#' Calculates the percentile for which to sum counts up to and scale by. Faster
#' version than available in cumNormStat. Deviates from methods described in Nature Methods by
#' making use of ro means for reference.
#' 
#' @param obj A MRexperiment object.
#' @param pFlag Plot the median difference quantiles.
#' @param rel Cutoff for the relative difference from one median difference
#' from the reference to the next.
#' @param ... Applicable if pFlag == TRUE. Additional plotting parameters.
#' @return Percentile for which to scale data
#' @seealso \code{\link{fitZig}} \code{\link{cumNorm}} \code{\link{cumNormStat}}
#' @examples
#' 
#' data(mouseData)
#' p = round(cumNormStatFast(mouseData,pFlag=FALSE),digits=2)
#' 
cumNormStatFast <-function(obj,pFlag = FALSE,rel=.1,...){
	if(class(obj)=="MRexperiment"){
		mat = MRcounts(obj,norm=FALSE,log=FALSE)
	} else {
		stop("Object needs to be a MRexperiment object.")
    }
	smat = lapply(1:ncol(mat), function(i) {
	    sort(mat[which(mat[, i]>0),i], decreasing = TRUE)
	})
	leng = max(sapply(smat,length))
	smat2 = array(NA,dim=c(leng,ncol(mat)))
	for(i in 1:ncol(mat)){
		smat2[leng:(leng-length(smat[[i]])+1),i] = smat[[i]]
	}

	rmat2 = sapply(1:ncol(smat2),function(i){
		quantile(smat2[,i],p=seq(0,1,length.out=nrow(smat2)),na.rm=TRUE)
	})
	smat2[is.na(smat2)] = 0
	ref1 = rowMeans(smat2)

	ncols = ncol(rmat2)
	diffr = sapply(1:ncols, function(i) {
		ref1 - rmat2[,i]
	})
	diffr1=rowMedians(abs(diffr))
	if(pFlag==TRUE){
		plot(abs(diff(diffr1))/diffr1[-1],type="h",...)
		abline(h=rel)
		axis(1,at=seq(0,length(diffr1),length.out=5),labels = seq(0,1,length.out=5))
	}
	x= which(abs(diff(diffr1))/diffr1[-1] > rel)[1]/length(diffr1)
	if(x<=0.50){
		warning("Low quantile estimate. Default value being used.")
		x = 0.50
	}
	obj@expSummary$cumNormStat = x;
	return(x)
}
