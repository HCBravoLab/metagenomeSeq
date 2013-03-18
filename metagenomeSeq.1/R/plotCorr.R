plotCorr <- function(obj,n,log=TRUE,norm=TRUE,fun=cor,...) {
    if (log == TRUE) {
        if (norm == TRUE) {
            mat = log2(cumNormMat(obj) + 1)
        }
        else {
            mat = log2(MRcounts(obj) + 1)
        }
    }
    else {
        if (norm == TRUE) {
            mat = cumNormMat(obj)
        }
        else {
            mat = MRcounts(obj)
        }
    }
    otusToKeep <- which(rowSums(mat) > 0)
    otuVars = rowSds(mat[otusToKeep, ])
    otuIndices = otusToKeep[order(otuVars, decreasing = TRUE)[1:n]]
    mat2 = mat[otuIndices, ]
    cc = as.matrix(fun(t(mat2)))
    hc = hclust(dist(mat2))
    otuOrder = hc$order
    cc = cc[otuOrder, otuOrder]
    heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
    heatmap.2(t(cc), col = heatmapCols, ...)
    invisible()
}

