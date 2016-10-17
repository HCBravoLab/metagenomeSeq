#' Load a count dataset associated with a study.
#' 
#' Load a matrix of OTUs in a tab delimited format
#' 
#' 
#' @aliases loadMeta metagenomicLoader
#' @param file Path and filename of the actual data file.
#' @param sep File delimiter.
#' @return A list with objects 'counts' and 'taxa'.
#' @seealso \code{\link{loadPhenoData}}
#' @examples
#' 
#' dataDirectory <- system.file("extdata", package="metagenomeSeq")
#' lung = loadMeta(file.path(dataDirectory,"CHK_NAME.otus.count.csv"))
#' 
loadMeta <- function(file,sep="\t")
{
	dat2 <- read.table(file,header=FALSE,sep=sep,nrows=1,stringsAsFactors=FALSE);
	subjects <- as.character(dat2[1,-1]);
	classes <-c("character",rep("numeric",length(subjects)));
	dat3 <- read.table(file,header=FALSE,skip=1,sep=sep,colClasses=classes,row.names=1);
	colnames(dat3) = subjects
	taxa<- rownames(dat3);
	obj <- list(counts=as.data.frame(dat3), taxa=as.data.frame(taxa))
	return(obj);
}
