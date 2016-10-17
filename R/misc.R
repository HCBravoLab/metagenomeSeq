#' Table of features unique to a group
#' 
#' Creates a table of features, their index, number of positive samples in a group,
#' and the number of reads in a group. Can threshold features by a minimum no. of reads
#' or no. of samples. 
#' 
#' @param obj Either a MRexperiment object or matrix.
#' @param cl A vector representing assigning samples to a group.
#' @param nsamples The minimum number of positive samples.
#' @param nreads The minimum number of raw reads.
#' @return Table of features unique to a group
#' @examples
#' data(mouseData)
#' head(uniqueFeatures(mouseData[1:100,],cl=pData(mouseData)[,3]))
#' 
uniqueFeatures<-function(obj,cl,nsamples=0,nreads=0){
	if (class(obj) == "MRexperiment") {
		mat = MRcounts(obj, norm = FALSE, log = FALSE)
	}
	else if (class(obj) == "matrix") {
		mat = obj
	}
	else {
		stop("Object needs to be either a MRexperiment object or matrix")
	}
	res = by(t(mat),cl,colSums)
	res = do.call("rbind",res)
	kreads = (colSums(res==0)>0)
	
	mat = mat>0
	resPos = by(t(mat),cl,colSums)
	resPos = do.call("rbind",resPos)
	ksamples = (colSums(resPos==0)>0)
	
	featureIndices = intersect(which(ksamples),which(kreads))
	numberReads = t(res[,featureIndices])
	colnames(numberReads) = paste("Reads in",colnames(numberReads))
	numberPosSamples = t(resPos[,featureIndices])
	colnames(numberPosSamples) = paste("Samp. in",colnames(numberPosSamples))
	featureIndices = featureIndices
	featureNames = rownames(mat[featureIndices,])
	
	df = cbind(featureIndices,numberPosSamples,numberReads)
	interesting = which(rowSums(numberReads)>=nreads & rowSums(numberPosSamples)>=nsamples)
	df[interesting,]
}

#' Function to make labels simpler
#' 
#' Beginning to transition to better axes for plots
#' 
#' @param x string for the x-axis
#' @param y string for the y-axis
#' @param norm is the data normalized?
#' @param log is the data logged?
#' @return vector of x,y labels
#' @examples
#' metagenomeSeq::makeLabels(norm=TRUE,log=TRUE) 
makeLabels<-function(x="samples",y="abundance",norm,log){
	yl = xl = ""
	if(log == TRUE){
		yl =  paste(yl,"Log2")
	}
	if(norm == TRUE){
		yl = paste(yl,"normalized")
	}
	yl = paste(yl,y)
	xl = paste(xl,x)
	return(c(xl,yl))
}
