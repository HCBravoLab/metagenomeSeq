#' Computes a slightly modified form of Metastats.
#' 
#' Wrapper to perform the permutation test on the t-statistic. This is the original
#' method employed by metastats (for non-sparse large samples). We include CSS normalization
#' though (optional) and log2 transform the data. In this method the null distribution is not assumed to be a t-dist.
#' 
#' 
#' @param obj A MRexperiment object with count data.
#' @param mod The model for the count distribution.
#' @param useCSSoffset Boolean, whether to include the default scaling
#' parameters in the model or not.
#' @param B Number of permutations.
#' @param coef The coefficient of interest.
#' @param sl The value to scale by (default=1000).
#'
#' @return Call made, fit object from lmFit, t-statistics and p-values for each feature.
#' @export
#' 
fitMeta <- function(obj,mod,useCSSoffset=TRUE,B=1000,coef=2,sl=1000){
    if(class(obj)=="MRexperiment"){
        mat = MRcounts(obj,norm=FALSE,log=FALSE)
        mat = log2(mat + 1)
    } else if(class(obj) == "matrix") {
        mat = obj
    } else {
        stop("Object needs to be either a MRexperiment object or matrix")
    }

	if(useCSSoffset==TRUE){
		if(any(is.na(normFactors(obj)))){
            stop("Calculate the normalization factors first!")
        }
		mmCount=cbind(mod,log2(normFactors(obj)/sl +1))}
	else{ 
       	mmCount=mod
   	}
    
    # fit of the data
	fitRes = limma::lmFit(mat,mmCount)	

    # The t-statistic
    tt <- fitRes$coef[,coef] / fitRes$stdev.unscaled[,coef] / fitRes$sigma

    perms = replicate(B,sample(mmCount[,coef]))
    mmCount1=mmCount[,-coef]

    nc = (ncol(mmCount1)+1)
    tobs<- sapply(1:B,function(i){
        # This code forces the covariate of interest to be a factor (might not apply)
        mmCountPerm = cbind(mmCount1,factor(perms[,i]))
        fit = limma::lmFit(mat,mmCountPerm)
        ttObs <- fit$coef[,nc] / fit$stdev.unscaled[,nc] / fit$sigma
        ttObs
    })
    p = rowMeans(abs(tobs)>=abs(tt))

	dat = list(call=match.call(),fit=fitRes,t = tt,p = p,type="perm")
	return(dat)
}
