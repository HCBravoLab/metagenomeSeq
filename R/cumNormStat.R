#' Cumulative sum scaling percentile selection
#' 
#' Calculates the percentile for which to sum counts up to and scale by.
#' cumNormStat might be deprecated one day. Deviates from methods in Nature Methods paper
#' by making use row means for generating reference.
#' 
#' @param obj A matrix or MRexperiment object.
#' @param qFlag Flag to either calculate the proper percentile using
#' R's step-wise quantile function or approximate function. 
#' @param pFlag Plot the relative difference of the median deviance from the reference.
#' @param rel Cutoff for the relative difference from one median difference
#' from the reference to the next
#' @param ... Applicable if pFlag == TRUE. Additional plotting parameters.
#' @return Percentile for which to scale data
#' @seealso \code{\link{fitZig}} \code{\link{cumNorm}} \code{\link{cumNormStatFast}}
#' @examples
#' 
#' data(mouseData)
#' p = round(cumNormStat(mouseData,pFlag=FALSE),digits=2)
#' 
cumNormStat <-
function(obj,qFlag = TRUE,pFlag = FALSE,rel=.1,...){
	mat = returnAppropriateObj(obj,FALSE,FALSE)
    if(any(colSums(mat)==0)) stop("Warning empty sample")
    
	smat = sapply(1:ncol(mat),function(i){sort(mat[,i],decreasing=FALSE)})
	ref  = rowMeans(smat);
	
	yy = mat;
	yy[yy==0]=NA;
	
	ncols = ncol(mat);
	refS = sort(ref);
    
	k = which(refS>0)[1]
	lo = (length(refS)-k+1)

	if(qFlag == TRUE){
		diffr = sapply(1:ncols,function(i){
		refS[k:length(refS)] - quantile(yy[,i],p=seq(0,1,length.out=lo),na.rm=TRUE)
	})
	}
	if(qFlag == FALSE){
       	diffr = sapply(1:ncols,function(i){
			refS[k:length(refS)] - approx(sort(yy[,i],decreasing=FALSE),n=lo)$y
		})
	}
	diffr2 = rowMedians(abs(diffr),na.rm=TRUE)
	if(pFlag ==TRUE){
        plot(abs(diff(diffr2[diffr2>0]))/diffr2[diffr2>0][-1],type="h",ylab="Relative difference for reference",xaxt="n",...)
		abline(h=rel)
		axis(1,at=seq(0,length(diffr2),length.out=5),labels = seq(0,1,length.out=5))
	}
	x = which(abs(diff(diffr2))/diffr2[-1]>rel)[1] / length(diffr2)
	if(x<=0.50){
		message("Default value being used.")
		x = 0.50
	}
	if(class("obj"=="MRexperiment")){
		obj@expSummary$cumNormStat = x;		
	}
	return(x)
}
