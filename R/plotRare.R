#' Plot of rarefaction effect
#' 
#' This function plots the number of observed features vs. the depth of coverage.
#' 
#' @param obj A MRexperiment object with count data or matrix.
#' @param cl Vector of classes for various samples.
#' @param ret True/False, return the number of features and the depth of coverage as a vector.
#' @param ... Additional plot arguments.
#' @return NA
#' @seealso \code{\link{plotOrd}}, \code{\link{plotMRheatmap}}, \code{\link{plotCorr}}, \code{\link{plotOTU}}, \code{\link{plotGenus}}
#' @examples
#' 
#' data(mouseData)
#' cl = factor(pData(mouseData)[,3])
#' res = plotRare(mouseData,cl=cl,ret=TRUE,pch=21,bg=cl)
#' tmp=lapply(levels(cl), function(lv) lm(res[,"ident"]~res[,"total"]-1, subset=cl==lv))
#' for(i in 1:length(levels(cl))){
#'    abline(tmp[[i]], col=i)
#' }
#' legend("topleft", c("Diet 1","Diet 2"), text.col=c(1,2),box.col=NA)
plotRare<-function(obj,cl=NULL,ret=FALSE,...){
    if(class(obj)=="MRexperiment"){
        mat = MRcounts(obj,norm=FALSE,log=FALSE)
        totalCounts = libSize(obj)
    } else if(class(obj) == "matrix") {
        mat = obj
        totalCounts=colSums(mat)
    } else {
        stop("Object needs to be either a MRexperiment object or matrix")
    }
    numFeatures=colSums(mat!=0)
    
    if(is.null(cl)){
        plot(totalCounts, numFeatures, xlab = "Depth of coverage", ylab = "Number of detected features",...)
    } else{
        plot(totalCounts, numFeatures, xlab = "Depth of coverage", ylab = "Number of detected features",col=factor(cl),...)
    }
    
    if(ret==TRUE){
        dat = cbind(totalCounts,numFeatures);
        colnames(dat) = c("total","ident")
        return(dat);
    }
    else{
        invisible()
    }
}




