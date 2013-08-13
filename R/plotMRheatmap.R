#' Basic heatmap plot function for normalized counts.
#' 
#' This function plots a heatmap of the "n" features with greatest variance
#' across rows.
#' 
#' 
#' @param obj A MRexperiment object with count data.
#' @param n The number of features to plot
#' @param trials A vector of clinical information for.
#' @param log Whether or not to log transform the counts.
#' @param norm Whether or not to normalize the counts.
#' @param ... Additional plot arguments.
#' @return NA
#' @seealso \code{\link{cumNormMat}}
#' @examples
#' 
#' data(mouseData)
#' trials = pData(mouseData)$diet
#' heatmapColColors=brewer.pal(12,"Set3")[as.integer(factor(trials))];
#' heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
#' plotMRheatmap(obj=mouseData,n=200,trials=trials,cexRow = 0.4,cexCol = 0.4,trace="none",
#'              col = heatmapCols,ColSideColors = heatmapColColors)
#' 
plotMRheatmap <- function(obj,n,trials,log=TRUE,norm=TRUE,...) {
  
  if(log==TRUE){
    if(norm==TRUE){
        mat = log2(cumNormMat(obj)+1);
    }else{
        mat = log2(MRcounts(obj)+1)
    }
  } else{
    if(norm==TRUE){
        mat = cumNormMat(obj)    
    }else{
        mat = MRcounts(obj)
    }
  }
  otusToKeep <- which(rowSums(mat)>0);
  otuVars=rowSds(mat[otusToKeep,]);
  otuIndices=otusToKeep[order(otuVars,decreasing=TRUE)[1:n]];
  mat2=mat[otuIndices,];
  heatmap.2(mat2,...);
  invisible()
}
