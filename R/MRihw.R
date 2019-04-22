#' MRihw runs IHW within a MRcoefs() call
#' 
#' Function used in MRcoefs() when "IHW" is set as the p value adjustment method
#' 
#' @rdname MRihw
#' @param obj Either a fitFeatureModelResults or fitZigResults object
#' @param ... other parameters
#' 
setGeneric("MRihw", function(obj, ...){standardGeneric("MRihw")})

#' MRihw runs IHW within a MRcoefs() call
#' 
#' Function used in MRcoefs() when "IHW" is set as the p value adjustment method
#' 
#' @rdname MRihw-fitFeatureModelResults
#' @param obj Either a fitFeatureModelResults or fitZigResults object
#' @param p a vector of pvalues extracted from obj
#' @param adjustMethod Value specifying which adjustment method and which covariate to use for IHW pvalue adjustment. 
#' For obj of class \code{\link{fitFeatureModelResults}}, options are "ihw-abundance" (median feature count per row) 
#' and "ihw-ubiquity" (number of non-zero features per row). For obj of class \code{\link{fitZigResults}}, 
#' options are "ihw-abundance" (weighted mean per feature) and "ihw-ubiquity" (number of non-zero features per row). 
#' @param alpha pvalue significance level specified for IHW call. Default is 0.1
#' 
setMethod("MRihw", signature = "fitFeatureModelResults", function(obj, p, adjustMethod, alpha){
  if (adjustMethod == "ihw-ubiquity") {
    # set covariate to be num of non-zero elements per row
    p <- obj@pvalues
    covariate <- rowSums(obj@counts != 0)
    ihwRes <- ihw(p, covariate, alpha)
    padj <- ihwRes@df$adj_pvalue 
  }
  if (adjustMethod == "ihw-abundance"){
    # use feature median count as covariate
    covariate <- rowMedians(obj@counts)
    ihwRes <- ihw(p, covariate, alpha)
    padj <- ihwRes@df$adj_pvalue
  }
  padj
})

#' MRihw runs IHW within a MRcoefs() call
#' 
#' Function used in MRcoefs() when "IHW" is set as the p value adjustment method
#' 
#' @rdname MRihw-fitZigResults
#' @param obj Either a fitFeatureModelResults or fitZigResults object
#' @param p a vector of pvalues extracted from obj
#' @param adjustMethod Value specifying which adjustment method and which covariate to use for IHW pvalue adjustment. 
#' For obj of class \code{\link{fitFeatureModelResults}}, options are "ihw-abundance" (median feature count per row) 
#' and "ihw-ubiquity" (number of non-zero features per row). For obj of class \code{\link{fitZigResults}}, 
#' options are "ihw-abundance" (weighted mean per feature) and "ihw-ubiquity" (number of non-zero features per row). 
#' @param alpha pvalue significance level specified for IHW call. Default is 0.1
#' 
setMethod("MRihw", signature = "fitZigResults", function(obj, p, adjustMethod, alpha){
  if (adjustMethod == "ihw-ubiquity"){
    #use number of non-zero features per row as the covariate in ihw() call
    covariate <- rowSums(obj@counts != 0)
    ihwRes <- ihw(p, covariate, alpha)
    padj <- ihwRes@df$adj_pvalue
  }
  if (adjustMethod == "ihw-abundance"){
    # use Amean as covariate
    covariate <- obj@eb$Amean
    ihwRes <- ihw(p, covariate, alpha)
    padj <- ihwRes@df$adj_pvalue
  }
  padj
})