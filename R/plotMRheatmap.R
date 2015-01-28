#' Basic heatmap plot function for normalized counts.
#' 
#' This function plots a heatmap of the 'n' features with greatest variance
#' across rows (or other statistic).
#' 
#' 
#' @param obj A MRexperiment object with count data.
#' @param n The number of features to plot. This chooses the 'n' features of greatest positive statistic.
#' @param norm Whether or not to normalize the counts - if MRexperiment object.
#' @param log Whether or not to log2 transform the counts - if MRexperiment object.
#' @param fun Function to select top 'n' features.
#' @param ... Additional plot arguments.
#' @return NA
#' @seealso \code{\link{cumNormMat}}
#' @examples
#' 
#' data(mouseData)
#' trials = pData(mouseData)$diet
#' heatmapColColors=brewer.pal(12,"Set3")[as.integer(factor(trials))];
#' heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
#' #### version using sd
#' plotMRheatmap(obj=mouseData,n=200,cexRow = 0.4,cexCol = 0.4,trace="none",
#'              col = heatmapCols,ColSideColors = heatmapColColors)
#' #### version using MAD
#' plotMRheatmap(obj=mouseData,n=50,fun=mad,cexRow = 0.4,cexCol = 0.4,trace="none",
#'              col = heatmapCols,ColSideColors = heatmapColColors) 
#'
plotMRheatmap <- function(obj,n,norm=TRUE,log=TRUE,fun=sd,...) {
  mat = returnAppropriateObj(obj,norm,log)
    
  otusToKeep = which(rowSums(mat)>0);
  otuStats   = apply(mat[otusToKeep,],1,fun);
  otuIndices = otusToKeep[order(otuStats,decreasing=TRUE)[1:n]];
  mat2=mat[otuIndices,];
  heatmap.2(mat2,...);
  invisible()
}
