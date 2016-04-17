#' Aggregates a MRexperiment object or counts matrix to by a factor.
#' 
#' Using the phenoData information in the MRexperiment, calling aggregateBySample on a
#' MRexperiment and a particular phenoData column (i.e. 'diet') will aggregate counts
#' using the aggfun function (default rowMeans). Possible aggfun alternatives
#' include rowMeans and rowMedians.
#' 
#' @param obj A MRexperiment object or count matrix.
#' @param fct phenoData column name from the MRexperiment object or if count matrix object a vector of labels.
#' @param aggfun Aggregation function.
#' @param out Either 'MRexperiment' or 'matrix'
#' @return An aggregated count matrix or MRexperiment object where the new pData is a vector of `fct` levels.
#' @aliases aggSamp
#' @rdname aggregateBySample
#' @export
#' @examples
#' 
#' data(mouseData)
#' aggregateBySample(mouseData[1:100,],fct="diet",aggfun=rowSums)
#' # not run
#' # aggregateBySample(mouseData,fct="diet",aggfun=matrixStats::rowMedians)
#' # aggSamp(mouseData,fct='diet',aggfun=rowMaxs)
#' 
aggregateBySample<-function(obj,fct,aggfun=rowMeans,out="MRexperiment"){
	if(class(obj)=="MRexperiment"){
		mat = MRcounts(obj)
		if(length(fct)==1) factors = as.character(pData(obj)[,fct])
		else factors = as.character(fct)
	} else {
		mat = obj
		factors = as.character(fct)
		if(length(factors)!=ncol(mat)) stop("If input is a count matrix, fct must be a vector of length = ncol(count matrix)")
	}
	if(!(out%in%c("MRexperiment","matrix"))){
		stop("The variable out must either be 'MRexperiment' or 'matrix'")
	}
	grps = split(seq_along(factors),factors)
	
	newMat = array(NA,dim=c(nrow(obj),length(grps)))
	for(i in seq_along(grps)){
		newMat[,i] = aggfun(mat[,grps[[i]],drop=FALSE])
	}
	colnames(newMat) = names(grps)
	rownames(newMat) = rownames(obj)
	if(out=='matrix') return(newMat)
	if(out=='MRexperiment'){
		pd = data.frame(names(grps))
		colnames(pd) = "phenoData"
		rownames(pd) = names(grps)
		pd = as(pd,"AnnotatedDataFrame")
		if(class(obj)=="MRexperiment"){
			fd = as(fData(obj),"AnnotatedDataFrame")
			newObj = newMRexperiment(newMat,featureData=fd,phenoData=pd)
		} else {
			newObj = newMRexperiment(newMat,phenoData=pd)
		}
		return(newObj)
	}
}
#' @rdname aggregateBySample
#' @export
aggSamp<-function(obj,fct,aggfun=rowMeans,out='MRexperiment'){
	aggregateBySample(obj,fct,aggfun=aggfun,out=out)
}
