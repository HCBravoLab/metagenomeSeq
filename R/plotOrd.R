#' Plot of either PCA or MDS coordinates for the distances of normalized or unnormalized counts.
#' 
#' This function plots the PCA / MDS coordinates for the "n" features of interest. Potentially uncovering batch
#' effects or feature relationships.
#' 
#' 
#' @param obj A MRexperiment object or count matrix.
#' @param tran Transpose the matrix.
#' @param comp Which components to display
#' @param usePCA TRUE/FALSE whether to use PCA  or MDS coordinates (TRUE is PCA).
#' @param useDist TRUE/FALSE whether to calculate distances.
#' @param dist.method If useDist==TRUE, what method to calculate distances.
#' @param log Whether or not to log2 the counts - if MRexperiment object.
#' @param norm Whether or not to normalize the counts - if MRexperiment object.
#' @param ret Whether or not to output the coordinates.
#' @param ntop Number of features to make use of in calculating your distances.
#' @param ... Additional plot arguments.
#' @return NA
#' @seealso \code{\link{cumNormMat}}
#' @examples
#' 
#' data(mouseData)
#' cl = pData(mouseData)[,3]
#' plotOrd(mouseData,tran=TRUE,useDist=TRUE,pch=21,bg=factor(cl),usePCA=FALSE)
#' 
plotOrd<-function(obj,tran=FALSE,comp=1:2,log=TRUE,norm=TRUE,usePCA=TRUE,useDist=FALSE,dist.method="euclidian",ret=FALSE,ntop=NULL,...){
    if(class(obj)=="MRexperiment"){
        mat = MRcounts(obj,norm=norm,log=log)
    } else if(class(obj) == "matrix") {
        mat = obj
    } else {
        stop("Object needs to be either a MRexperiment object or matrix")
    }
    if(useDist==FALSE & usePCA==FALSE) stop("Classical MDS requires distances")
    if(is.null(ntop)) ntop = min(nrow(mat),1000)
    if(length(comp)>2) stop("Can't display more than two components")

    otusToKeep <- which(rowSums(mat)>0)
    otuVars<-rowSds(mat[otusToKeep,])
    otuIndices<-otusToKeep[order(otuVars,decreasing=TRUE)[seq_len(ntop)]]
    mat <- mat[otuIndices,]

    if(tran==TRUE){
        mat = t(mat)
    }
    if(useDist==TRUE){
        d <- dist(mat,method=dist.method)
    } else{ d = mat }
    
    if(usePCA==FALSE){
	   ord = cmdscale(d,k = max(comp))
    } else{
       ord = prcomp(d)$x
    }
    main = paste(ifelse(usePCA,"PCA","MDS"),ifelse(useDist,"distances of features","features"),sep=" on ")

	plot(ord[,comp],main=main,...)
    if(ret==TRUE) return(ord[,comp])
    else invisible()
}
