#' @name aggregateByTaxonomy
#' @title Aggregates a MRexperiment object by a particular taxonomic level.
#' 
#' Using the featureData information in the MRexperiment, calling aggregateByTaxonomy on a
#' MRexperiment and a particular featureData column (i.e. 'genus') will aggregate counts
#' to the desired level by with the aggfun function (default colSums). Possible aggfun alternatives
#' include colMeans and colMedians.
#' 
#' @param obj A MRexperiment object.
#' @param lvl featureData column name from the MRexperiment object.
#' @param norm Whether to aggregate normalized counts or not.
#' @param aggfun Aggregation function.
#' @return An aggregated count matrix.
#' @aliases aggTax
#' @rdname aggregateByTaxonomy
#' @export
#' @examples
#' 
#' # not run
#' # aggregateByTaxonomy(mouseData,lvl="genus",norm=TRUE,aggfun=colMedians)
#' 
aggregateByTaxonomy<-function(obj,lvl,norm=TRUE,aggfun = colSums){
    if(class(obj)=="MRexperiment"){
        mat = MRcounts(obj,norm=norm,log=FALSE)
    } else {
        stop("Object needs to be a MRexperiment object. If it's a matrix, see aggregateM.")
    }
    
	levels = fData(obj)[,lvl]
	grps = split(seq_along(levels),levels)
	
	newMat = array(NA,dim=c(length(grps),ncol(obj)))
	for(i in 1:length(grps)){
		newMat[i,] = aggfun(mat[grps[[i]],,drop=FALSE])
	}
	rownames(newMat) = names(grps)
	colnames(newMat) = sampleNames(obj)
	newMat
}
#' @rdname aggregateByTaxonomy
#' @export
aggTax<-function(obj,lvl,norm=TRUE,aggfun = colSums){
    aggregateByTaxonomy(obj,lvl,norm=TRUE,aggfun = colSums)
}