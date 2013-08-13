#' Various statistics of the count data.
#' 
#' A matrix of values for each sample. The matrix consists of sample ids, the
#' sample scaling factor, quantile value, and the number of number of features.
#' 
#' 
#' @param obj A MRexperiment object with count data.
#' @param p Quantile value to calculate the scaling factor and quantiles for
#' the various samples.
#' @param output Output file name.
#' @return None.
#' @seealso \code{\link{cumNorm}} \code{\link{quantile}}
#' @examples
#' 
#' # see vignette
#' 
exportStats <-
function(obj,p= cumNormStat(obj),output="~/Desktop/res.stats.tsv"){

			xx=MRcounts(obj)
			xx[xx==0]=NA
			qs=matrixStats::colQuantiles(xx,p=p,na.rm=TRUE)
			s95 = colSums(xx<=qs,na.rm=TRUE)
			xx[xx>0] = 1;
			xx[is.na(xx)]=0
			
			newMat <- array(NA,dim=c(4,ncol(xx)+1));
			newMat[1,1] = "Subject"
			newMat[2,1] = "Scaling factor"
			newMat[3,1] = "Quantile value"
			newMat[4,1] = "Number of features"

			newMat[1,2:ncol(newMat)]<-sampleNames(obj);
			newMat[2,2:ncol(newMat)]<-s95;
			newMat[3,2:ncol(newMat)]<-qs;
			newMat[4,2:ncol(newMat)]<-colSums(xx)

			write((newMat),file = output,sep = "\t",ncolumns = 4);
}
