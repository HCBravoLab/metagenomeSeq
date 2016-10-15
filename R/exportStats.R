#' Various statistics of the count data.
#' 
#' A matrix of values for each sample. The matrix consists of sample ids, the
#' sample scaling factor, quantile value, the number identified features, and library size (depth of coverage).
#' 
#' 
#' @param obj A MRexperiment object with count data.
#' @param p Quantile value to calculate the scaling factor and quantiles for
#' the various samples.
#' @param file Output file name.
#' @return None.
#' @seealso \code{\link{cumNorm}} \code{\link{quantile}}
#' @examples
#' 
#' data(lungData)
#' dataDirectory <- system.file("extdata", package="metagenomeSeq")
#' exportStats(lungData[,1:5],file=file.path(dataDirectory,"tmp.tsv"))
#' head(read.csv(file=file.path(dataDirectory,"tmp.tsv"),sep="\t"))
#' 
exportStats <-function(obj,p= cumNormStat(obj),file="~/Desktop/res.stats.tsv"){
	xx=MRcounts(obj)
	xx[xx==0]=NA
	qs=colQuantiles(xx,probs=p,na.rm=TRUE)

	xx[xx>0] = 1;
	xx[is.na(xx)]=0
			
	newMat <- array(NA,dim=c(5,ncol(xx)+1));
	newMat[1,1] = "Subject"
	newMat[2,1] = "Scaling factor"
	newMat[3,1] = "Quantile value"
	newMat[4,1] = "Number of identified features"
	newMat[5,1] = "Library size"

	newMat[1,2:ncol(newMat)]<-sampleNames(obj);
	newMat[2,2:ncol(newMat)]<-unlist(normFactors(obj));
	newMat[3,2:ncol(newMat)]<-qs;
	newMat[4,2:ncol(newMat)]<-colSums(xx);
	newMat[5,2:ncol(newMat)]<-unlist(libSize(obj));

	write((newMat),file = file,sep = "\t",ncolumns = 5);
}
