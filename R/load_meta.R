#' Load a count dataset associated with a study.
#' 
#' Load a matrix of OTUs in a tab delimited format
#' 
#' 
#' @aliases load_meta metagenomicLoader
#' @param file Path and filename of the actual data file.
#' @param sep File delimiter.
#' @return An object of count data.
#' @seealso \code{\link{load_phenoData}}
#' @examples
#' 
#' dataDirectory <- system.file("extdata", package="metagenomeSeq")
#' lung = load_meta(file.path(dataDirectory,"CHK_NAME.otus.count.csv"))
#' 
load_meta <- function(file,sep="\t")
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
