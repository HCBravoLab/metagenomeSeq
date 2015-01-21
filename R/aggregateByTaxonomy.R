#' @name aggregateByTaxonomy
#' @title Aggregates a MRexperiment object or counts matrix to a particular level.
#' 
#' @details Using the featureData information in the MRexperiment, calling aggregateByTaxonomy on a
#' MRexperiment and a particular featureData column (i.e. 'genus') will aggregate counts
#' to the desired level using the aggfun function (default colSums). Possible aggfun alternatives
#' include colMeans and colMedians.
#' 
#' @param obj A MRexperiment object or count matrix.
#' @param lvl featureData column name from the MRexperiment object or if count matrix object a vector of labels.
#' @param alternate Use the rowname for undefined OTUs instead of aggregating to "no_match".
#' @param norm Whether to aggregate normalized counts or not.
#' @param log Whether or not to log2 transform the counts - if MRexperiment object.
#' @param aggfun Aggregation function.
#' @param sl scaling value, default is 1000.
#' @param out Either 'MRexperiment' or 'matrix'
#' @return An aggregated count matrix.
#' @aliases aggTax
#' @rdname aggregateByTaxonomy
#' @export
#' @examples
#' 
#' data(mouseData)
#' aggregateByTaxonomy(mouseData[1:100,],lvl="class",norm=TRUE,aggfun=colSums)
#' # not run
#' # aggregateByTaxonomy(mouseData,lvl="class",norm=TRUE,aggfun=colMedians)
#' # aggTax(mouseData,lvl='phylum',norm=FALSE,aggfun=colSums)
#' 
aggregateByTaxonomy<-function(obj,lvl,alternate=FALSE,norm=FALSE,log=FALSE,aggfun = colSums,sl=1000,out="MRexperiment"){
	if(class(obj)=="MRexperiment"){
		mat = MRcounts(obj,norm=norm,log=log,sl=sl)
		if(length(lvl)==1) levels = as.character(fData(obj)[,lvl])
		else levels = as.character(lvl)
	} else {
		mat = obj
		levels = as.character(lvl)
		if(length(levels)!=nrow(mat)) stop("If input is a count matrix, lvl must be a vector of length = nrow(count matrix)")
	}
	if(!(out%in%c("MRexperiment","matrix"))){
		stop("The variable out must either be 'MRexperiment' or 'matrix'")
	}
	
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
	for(i in seq_along(grps)){
		newMat[i,] = aggfun(mat[grps[[i]],,drop=FALSE])
	}
	rownames(newMat) = names(grps)
	colnames(newMat) = colnames(obj)
	if(out=='matrix') return(newMat)
	if(out=='MRexperiment'){
		taxa = data.frame(names(grps))
		colnames(taxa) = "Taxa"
		rownames(taxa) = names(grps)
		taxa = as(taxa,"AnnotatedDataFrame")
		if(class(obj)=="MRexperiment"){
			pd = as(pData(obj),"AnnotatedDataFrame")
			newObj = newMRexperiment(newMat,featureData=taxa,phenoData=pd)
		} else {
			newObj = newMRexperiment(newMat,featureData=taxa)
		}
		return(newObj)
	}
}
#' @rdname aggregateByTaxonomy
#' @export
aggTax<-function(obj,lvl,alternate=FALSE,norm=FALSE,log=FALSE,aggfun = colSums,sl=1000,out='MRexperiment'){
	aggregateByTaxonomy(obj,lvl,alternate=alternate,norm=norm,log=log,aggfun = aggfun,sl=sl,out=out)
}
