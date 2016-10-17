#' Extract the essentials of an MRexperiment.
#'
#' @param obj MRexperiment-class object.
#'
#' @return \itemize{A list containing:
#' \item     counts        : Count data
#' \item     librarySize         : The column sums / library size / sequencing depth 
#' \item     normFactors          : The normalization scaling factors
#' \item     pheno  : phenotype table
#' \item     feat  : feature table
#' }
#' 
#' @examples
#'
#' data(mouseData)
#' head(metagenomeSeq:::extractMR(mouseData))
#'
extractMR<-function(obj){
	mat = MRcounts(obj)
	ls  = as.vector(libSize(obj))
	norm= as.vector(normFactors(obj))
	pd  = pData(obj)
	fd  = fData(obj)
	dat = list(counts=mat,librarySize=ls,normFactors=norm,pheno=pd,feat=fd)
	return(dat)
}

#' Merge two tables
#'
#' @param x Table 1.
#' @param y Table 2.
#'
#' @return Merged table
#'
mergeTable<-function(x,y){
	rows = union(rownames(x),rownames(y))
	cols = union(colnames(x),colnames(y))
	fullmat = array(NA,dim=c(length(rows),length(cols)))
	rownames(fullmat) = rows
	colnames(fullmat) = cols
	fullmat[rownames(x),colnames(x)] = as.matrix(x)
	fullmat[rownames(y),colnames(y)] = as.matrix(y)
	fullmat
}

#' Merge two MRexperiment objects together
#'
#' This function will take two MRexperiment objects and merge them together finding common
#' OTUs. If there are OTUs not found in one of the two MRexperiments then a message will
#' announce this and values will be coerced to zero for the second table.
#'
#' @param x MRexperiment-class object 1.
#' @param y MRexperiment-class object 2.
#'
#' @return Merged MRexperiment-class object.
#' @export
#'
#' @examples
#' data(mouseData)
#' newobj = mergeMRexperiments(mouseData,mouseData)
#' newobj
#' 
#' # let me know if people are interested in an option to merge by keys instead of row names.
#' data(lungData)
#' newobj = mergeMRexperiments(mouseData,lungData)
#' newobj
mergeMRexperiments<-function(x,y){
	xdat = extractMR(x)
	ydat = extractMR(y)
	xmat = xdat$counts; ymat = ydat$counts
	cnames = union(colnames(xmat),colnames(ymat))
	if(length(cnames)!=(ncol(x)+ncol(y))){
		message("MRexperiment 1 and 2 share sample ids; adding labels to sample ids.")
		newXnames = paste(colnames(xmat),"x",sep=".")
		newYnames = paste(colnames(ymat),"y",sep=".")
		cnames = union(newXnames,newYnames)
		colnames(xdat$counts) = 
			rownames(xdat$pheno) = 
			names(xdat$normFactors) = 
			names(xdat$librarySize) = 
			newXnames
		colnames(ydat$counts) = 
			rownames(ydat$pheno) = 
			names(ydat$normFactors) = 
			names(ydat$librarySize) = 
			newYnames
	}
	counts = mergeTable(xdat$counts,ydat$counts)
	pheno  = as.data.frame(mergeTable(xdat$pheno,ydat$pheno))
	feat   = as.data.frame(mergeTable(xdat$feat,ydat$feat))
	librarySize = c(xdat$librarySize,ydat$librarySize)
	normFactors = c(xdat$normFactors,ydat$normFactors)
	if(any(is.na(counts))){
		message("There were OTUs not shared between objects. Coercing values to 0.")
		counts[is.na(counts)] = 0
	}
	obj = newMRexperiment(counts=counts,
					normFactors=normFactors,
					libSize=librarySize,
					phenoData = AnnotatedDataFrame(pheno),
					featureData=AnnotatedDataFrame(feat))
	return(obj)
}
