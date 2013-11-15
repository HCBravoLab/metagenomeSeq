#' Basic heatmap plot function for normalized counts.
#' 
#' This function plots a heatmap of the "n" features with greatest variance
#' across rows.
#' 
#' 
#' @param obj A MRexperiment object with count data.
#' @param n The number of features to plot. This chooses the "n" features with greatest variance.
#' @param log Whether or not to log2 transform the counts - if MRexperiment object.
#' @param norm Whether or not to normalize the counts - if MRexperiment object.
#' @param ... Additional plot arguments.
#' @return NA
#' @seealso \code{\link{cumNormMat}}
#' @examples
#' 
#' data(mouseData)
#' trials = pData(mouseData)$diet
#' heatmapColColors=brewer.pal(12,"Set3")[as.integer(factor(trials))];
#' heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
#' plotMRheatmap(obj=mouseData,n=200,cexRow = 0.4,cexCol = 0.4,trace="none",
#'              col = heatmapCols,ColSideColors = heatmapColColors)
#' 
plotMRheatmap <- function(obj,n,log=TRUE,norm=TRUE,...) {
    if(class(obj)=="MRexperiment"){
        mat = MRcounts(obj,norm=norm,log=log)
    } else if(class(obj) == "matrix") {
        mat = obj
    } else {
       stop("Object needs to be either a MRexperiment object or matrix")
    }
    
  otusToKeep <- which(rowSums(mat)>0);
  otuVars=rowSds(mat[otusToKeep,]);
  otuIndices=otusToKeep[order(otuVars,decreasing=TRUE)[1:n]];
  mat2=mat[otuIndices,];
  heatmap.2(mat2,...);
  invisible()
}
