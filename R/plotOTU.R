#' Basic plot function of the raw or normalized data.
#' 
#' This function plots the abundance of a particular OTU by class. The function
#' uses the estimated posterior probabilities to make technical zeros
#' transparent.
#' 
#' 
#' @param obj A MRexperiment object with count data.
#' @param otu The row number/OTU to plot.
#' @param classIndex A list of the samples in their respective groups.
#' @param log Whether or not to log2 transform the counts - if MRexperiment object.
#' @param norm Whether or not to normalize the counts - if MRexperiment object.
#' @param jitter.factor Factor value for jitter.
#' @param pch Standard pch value for the plot command.
#' @param labs Whether to include group labels or not. (TRUE/FALSE)
#' @param xlab xlabel for the plot.
#' @param ylab ylabel for the plot.
#' @param jitter Boolean to jitter the count data or not.
#' @param ... Additional plot arguments.
#' @return Plotted values
#' @seealso \code{\link{cumNorm}}
#' @examples
#' 
#' data(mouseData)
#' classIndex=list(controls=which(pData(mouseData)$diet=="BK"))
#' classIndex$cases=which(pData(mouseData)$diet=="Western")
#' # you can specify whether or not to normalize, and to what level
#' plotOTU(mouseData,otu=9083,classIndex,norm=FALSE,main="9083 feature abundances")
#' 
plotOTU <-
function(obj,otu,classIndex,log=TRUE,norm=TRUE,jitter.factor=1,pch=21,labs=TRUE,xlab=NULL,ylab=NULL,jitter=TRUE,...){
    mat = returnAppropriateObj(obj,norm,log)
    
	l=lapply(classIndex, function(j){
        mat[otu,j]
        })

	z = posteriorProbs(obj)
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
    
	if(jitter) x=jitter(x,jitter.factor)
    
    if(is.null(ylab)){ylab="Normalized log(cpt)"}
    if(is.null(xlab)){xlab="Groups of comparison"}

    plot(x,y,col=col,pch=pch,bg=col,xlab=xlab,ylab=ylab,xaxt="n",...)
    if(labs==TRUE){
        gp = names(classIndex)
        axis(1,at=seq(1:length(gp)),gp)
    }
	invisible(list(x=x,y=y))
}
