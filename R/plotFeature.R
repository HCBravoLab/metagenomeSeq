#' Basic plot function of the raw or normalized data.
#' 
#' This function plots the abundance of a particular OTU by class. The function
#' is the typical manhattan plot of the abundances.
#'
#' 
#' @param obj A MRexperiment object with count data.
#' @param otuIndex The row to plot
#' @param classIndex A list of the samples in their respective groups.
#' @param col A vector to color samples by.
#' @param sort Boolean, sort or not.
#' @param sortby Default is sort by library size, alternative vector for sorting
#' @param norm Whether or not to normalize the counts - if MRexperiment object.
#' @param log Whether or not to log2 transform the counts - if MRexperiment object.
#' @param sl Scaling factor - if MRexperiment and norm=TRUE.
#' @param ... Additional plot arguments.
#' @return counts and classindex
#' @seealso \code{\link{cumNorm}}
#' @examples
#' 
#' data(mouseData)
#' classIndex=list(Western=which(pData(mouseData)$diet=="Western"))
#' classIndex$BK=which(pData(mouseData)$diet=="BK")
#' otuIndex = 8770
#'
#' par(mfrow=c(2,1))
#' dates = pData(mouseData)$date
#' plotFeature(mouseData,norm=FALSE,log=FALSE,otuIndex,classIndex,
#' col=dates,sortby=dates,ylab="Raw reads")
#'
plotFeature<-function(obj,otuIndex,classIndex,col="black",sort=TRUE,sortby=NULL,norm=TRUE,log=TRUE,sl=1000,...){
    mat = returnAppropriateObj(obj,norm,log,sl)
    fmat = mat[otuIndex,]
    ylmin = min(fmat)
    ylmax = max(fmat)

    nplots = length(classIndex)
    nms = names(classIndex)

    counts = lapply(classIndex,function(i){
        fmat[i]
        })
    if(sort==TRUE){
        if(is.null(sortby)){
            ord  = lapply(classIndex,function(i){
                order(colSums(mat[,i]))
                })
        } else{
            ord  = lapply(classIndex,function(i){
                order(sortby[i])
                })
        }
    } else {
        ord = lapply(classIndex,function(i){
            1:length(i)
            })
    }
    if(length(col)>1){
        col = as.integer(factor(col))
        col4groups = lapply(1:length(classIndex),function(i){
                cindex = classIndex[[i]]
                oindex = ord[[i]]
                col[cindex[oindex]]
        })
    }
    for(i in 1:nplots){
        vals = counts[[i]][ord[[i]]]
        if(exists("col4groups")) colors  = col4groups[[i]]
        else colors = col
        plot(vals,xlab=nms[i],type="h",col=colors,ylim=c(ylmin,ylmax),...)
    }
    invisible(cbind(counts,ord))
}
