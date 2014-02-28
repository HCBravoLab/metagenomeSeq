#' @name aggregateByTaxonomy
#' @title Aggregates a MRexperiment object by a particular taxonomic level.
#' 
#' Using the featureData information in the MRexperiment, calling aggregateByTaxonomy on a
#' MRexperiment and a particular featureData column (i.e. 'genus') will aggregate counts
#' to the desired level using the aggfun function (default colSums). Possible aggfun alternatives
#' include colMeans and colMedians.
#' 
#' @param obj A MRexperiment object.
#' @param lvl featureData column name from the MRexperiment object.
#' @param alternate Use the rowname for undefined OTUs instead of aggregating to "no_match".
#' @param norm Whether to aggregate normalized counts or not.
#' @param log Whether or not to log2 transform the counts - if MRexperiment object.
#' @param aggfun Aggregation function.
#' @param sl scaling value, default is 1000.
#' @return An aggregated count matrix.
#' @aliases aggTax
#' @rdname aggregateByTaxonomy
#' @export
#' @examples
#' 
#' # not run
#' # aggregateByTaxonomy(mouseData,lvl="genus",norm=TRUE,aggfun=colMedians)
#' 
aggregateByTaxonomy<-function(obj,lvl,alternate=FALSE,norm=TRUE,log=FALSE,aggfun = colSums,sl=1000){
	if(class(obj)=="MRexperiment"){
		mat = MRcounts(obj,norm=norm,log=log,sl=sl)
    } else {
		stop("Object needs to be a MRexperiment object. If it's a matrix, see aggregateM.")
    }
    
	levels = as.character(fData(obj)[,lvl])
    nafeatures = is.na(levels)
    if(length(nafeatures)>0){
        if(alternate==FALSE){
            levels[nafeatures] = "no_match"
        } else {
            levels[nafeatures] = paste("OTU_",rownames(obj)[nafeatures],sep="")
        }
    }
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
aggTax<-function(obj,lvl,alternate=FALSE,norm=TRUE,log=FALSE,aggfun = colSums,sl=1000){
    aggregateByTaxonomy(obj,lvl,alternate=alternate,norm=norm,log=log,aggfun = aggfun,sl=sl)
}