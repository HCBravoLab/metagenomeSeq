#' Basic plot function of the raw or normalized data.
#' 
#' This function plots the abundance of a particular OTU by class. The function
#' uses the estimated posterior probabilities to make technical zeros
#' transparent.
#' 
#' 
#' @aliases genusPlot plotGenus
#' @param obj An MRexperiment object with count data.
#' @param otuIndex A list of the otus with the same annotation.
#' @param classIndex A list of the samples in their respective groups.
#' @param log Whether or not to log transform the counts - if MRexperiment object.
#' @param norm Whether or not to normalize the counts - if MRexperiment object.
#' @param no Which of the otuIndex to plot.
#' @param jitter.factor Factor value for jitter
#' @param pch Standard pch value for the plot command.
#' @param labs Whether to include group labels or not. (TRUE/FALSE)
#' @param xlab xlabel for the plot.
#' @param ylab ylabel for the plot.
#' @param jitter Boolean to jitter the count data or not.
#' @param ret Boolean to return the observed data that would have been plotted.
#' @param ... Additional plot arguments.
#' @return NA
#' @seealso \code{\link{cumNorm}}
#' @examples
#' 
#' data(mouseData)
#' classIndex=list(controls=which(pData(mouseData)$diet=="BK"))
#' classIndex$cases=which(pData(mouseData)$diet=="Western")
#' otuIndex = grep("Strep",fData(mouseData)$fdata)
#' otuIndex=otuIndex[order(rowSums(MRcounts(mouseData)[otuIndex,]),decreasing=TRUE)]
#' plotGenus(mouseData,otuIndex,classIndex,no=1:2,xaxt="n",norm=FALSE,ylab="Strep normalized log(cpt)")
#' 
plotGenus <-
function(obj,otuIndex,classIndex,log=TRUE,norm=TRUE,no=1:length(otuIndex),labs=TRUE,xlab=NULL,ylab=NULL,jitter=TRUE,jitter.factor=1,pch=21,ret=FALSE,...){

    if(class(obj)=="MRexperiment"){
        mat = MRcounts(obj,norm=norm,log=log)
    } else if(class(obj) == "matrix") {
        mat = obj
    } else {
       stop("Object needs to be either a MRexperiment object or matrix")
    }

	l=lapply(otuIndex[no], function(i) lapply(classIndex, function(j) {
        mat[i,j]
        }))

	l=unlist(l,recursive=FALSE)
	if(!is.list(l)) stop("l must be a list\n")
	y=unlist(l)
	x=rep(seq(along=l),sapply(l,length))

	z = posterior.probs(obj)
    #if(!is.null(z)){
    #    z = 1-z;
    #    lz=lapply(classIndex,function(j){(z[otuIndex[no],j])})
    #    z = unlist(lz)
    #    blackCol=t(col2rgb("black"))
    #    col=rgb(blackCol,alpha=z)
    #} else {
        blackCol=t(col2rgb("black"))
        col=rgb(blackCol)
    #}
    if(jitter) x=jitter(x,jitter.factor)
    
    if(is.null(ylab)){ylab="Normalized log(cpt)"}
    if(is.null(xlab)){xlab="Groups of comparison"}

    plot(x,y,col=col,pch=pch,xlab=xlab,ylab=ylab,xaxt="n",...)
    if(labs==TRUE){
        gp = rep(names(classIndex),length(no))
        axis(1,at=seq(1:length(gp)),gp)
    }

	if(ret) list(x=x,y=y)
}
