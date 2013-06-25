#' Cumulative normalization statistic.
#'
#' @param obj An eSet object.
#' @param pFlag Whether or not to plot the reference.
#' @param rel Relative difference of rel percent.
#' @return P-value for which to cumulative normalize.
#'
#' @name cumNormStat
#' @seealso \code{\link{fitZig}} \code{\link{cumNorm}}
cumNormStat <-
function(obj,pFlag = FALSE,rel=.1,qFlag = TRUE, ...){
    
	mat = MRcounts(obj);
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
            refS[k:length(refS)] - approx(yy[,i],n=lo)$y
        })
    }
	diffr2 = matrixStats::rowMedians(abs(diffr),na.rm=TRUE)
	if(pFlag ==TRUE){
        plot(abs(diff(diffr2[diffr2>0]))/diffr2[diffr2>0][-1],type="h",ylab="Relative difference for reference",xaxt="n",...)
		abline(h=rel)
        axis(1,at=seq(0,length(diffr2),length.out=5),labels = seq(0,1,length.out=5))
	}
    x = which(abs(diff(diffr2))/diffr2[-1]>rel)[1] / length(diffr2)
    obj@expSummary$cumNormStat = x;
	return(x)
}