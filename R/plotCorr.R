#' Basic correlation plot function for normalized or unnormalized counts.
#' 
#' This function plots a heatmap of the "n" features with greatest variance
#' across rows.
#' 
#' 
#' @param obj A MRexperiment object with count data.
#' @param n The number of features to plot. This chooses the "n" features with greatest variance.
#' @param norm Whether or not to normalize the counts - if MRexperiment object.
#' @param log Whether or not to log2 transform the counts - if MRexperiment object.
#' @param fun Function to calculate pair-wise relationships. Default is pearson
#' correlation
#' @param ... Additional plot arguments.
#' @return plotted correlation matrix
#' @seealso \code{\link{cumNormMat}}
#' @examples
#' 
#' data(mouseData)
#' plotCorr(obj=mouseData,n=200,cexRow = 0.4,cexCol = 0.4,trace="none",dendrogram="none",
#'          col = colorRampPalette(brewer.pal(9, "RdBu"))(50))
#' 
plotCorr <- function(obj,n,norm=TRUE,log=TRUE,fun=cor,...) {
    mat = returnAppropriateObj(obj,norm,log)
    otusToKeep <- which(rowSums(mat) > 0)
    otuVars = rowSds(mat[otusToKeep, ])
    otuIndices = otusToKeep[order(otuVars, decreasing = TRUE)[1:n]]
    mat2 = mat[otuIndices, ]
    cc = as.matrix(fun(t(mat2)))
    hc = hclust(dist(mat2))
    otuOrder = hc$order
    cc = cc[otuOrder, otuOrder]
    heatmap.2(t(cc),...)
    invisible(t(cc))
}

