#' Basic plot function of the raw or normalized data.
#'
#' This function plots the abundance of a particular OTU by class. The function uses
#' the estimated posterior probabilities to make technical zeros transparent. 
#'
#' @param obj An eSet object with count data.
#' @param otuIndex A list of the otus with the same annotation.
#' @param classIndex A list of the samples in their respective groups.
#' @param norm Whether or not to normalize the counts.
#' @param normp The value at which to scale the counts by and then log.
#' @param no Which of the otuIndex to plot.
#' @param factor Factor value for jitter
#' @param pch Standard pch value for the plot command.
#' @param jitter Boolean to jitter the count data or not.
#' @param ret Boolean to return the observed data that would have been plotted.
#' @param ... Additional plot arguments.
#' @return NA
#' @note \code{\link{detect}} makes use of settings.
#'
#' @name plotGenus
#' @aliases genusPlot
#' @seealso \code{\link{cumNorm}}
#' @examples 
#' classIndex=list(controls=which(type=="Control"))
#' classIndex$cases=which(type=="Case")
#' otuIndex = which(taxa == "E-coli")
#' plotGenus(obj,otu=12,classIndex,xlab="OTU log-normalized counts")

plotGenus <-
function(obj,otuIndex,classIndex,norm=TRUE,no=1:length(otuIndex),labs=TRUE,xlab=NULL,ylab=NULL,jitter=TRUE,factor=1,pch=21,ret=FALSE,...){

	l=lapply(otuIndex[no], function(i) lapply(classIndex, function(j) {
        if(norm==FALSE){ x=log2(MRcounts(obj)[i,j]+1) }
        else if(norm==TRUE){  x=log2(MRcounts(obj,norm=TRUE)+1)[i,j]; }
        x 
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
    if(jitter) x=jitter(x,factor)
    
    if(is.null(ylab)){ylab="Normalized log(cpt)"}
    if(is.null(xlab)){xlab="Groups of comparison"}

    plot(x,y,col=col,pch=pch,xlab=xlab,ylab=ylab,xaxt="n",...)
    if(labs==TRUE){
        gp = rep(names(classIndex),length(no))
        axis(1,at=seq(1:length(gp)),gp)
    }

	if(ret) list(x=x,y=y)
}