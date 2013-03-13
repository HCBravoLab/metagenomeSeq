#' Basic plot function of the raw or normalized data.
#'
#' This function plots the abundance of a particular OTU by class. The function uses
#' the estimated posterior probabilities to make technical zeros transparent. 
#'
#' @param obj An eSet object with count data.
#' @param otu The row number/OTU to plot.
#' @param classIndex A list of the samples in their respective groups.
#' @param norm Whether or not to normalize the counts.
#' @param normp The value at which to scale the counts by and then log.
#' @param factor Factor value for jitter.
#' @param pch Standard pch value for the plot command.
#' @param jitter Boolean to jitter the count data or not.
#' @param ret Boolean to return the observed data that would have been plotted.
#' @param ... Additional plot arguments.
#' @return NA
#'
#' @name plotOTU
#' @aliases otuplot
#' @seealso \code{\link{cumNorm}}
#' @examples 
#' classIndex=list(controls=which(type=="Control"))
#' classIndex$cases=which(type=="Case")
#' plotOTU(obj,otu=12,classIndex,xlab="OTU log-normalized counts")

plotOTU <-
function(obj,otu,classIndex,norm=TRUE,factor=1,pch=21,jitter=TRUE,ret=FALSE,...){

	l=lapply(classIndex, function(j){
        if(norm==FALSE){ 
            log2(MRcounts(obj)[otu,j]+1) 
        }
        else if(norm==TRUE){
            if(any(is.na(normFactors(obj)))){
                log2(cumNormMat(obj)[otu,j]+1)
            } else{ 
                log2(MRcounts(obj,norm=TRUE)[otu,j]+1)
            }
        }
        })

	z = posterior.probs(obj)
    y=unlist(l)
    x=rep(seq(along=l),sapply(l,length))

    if(!is.null(z)){
        z = 1-z;
        lz=lapply(classIndex,function(j){(z[otu,j])})
        z = unlist(lz)
        blackCol=t(col2rgb("black"))
        col=rgb(blackCol,alpha=z)
    } else {
        blackCol=t(col2rgb("black"))
        col=rgb(blackCol)
    }
    
	if(jitter) x=jitter(x,factor)

	plot(x,y,col=col,pch=pch,bg=col,...)
	if (ret)
		list(x=x,y=y)
}
