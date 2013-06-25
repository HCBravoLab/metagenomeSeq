#' Load a count dataset associated with a study.
#'
#' Load a matrix of OTUs in a tab delimited format
#'
#' @param file Path and filename of the actual data file.
#' @return An object of count data.
#'
#' @name load_meta
#' @aliases metagenomicLoader
#' @seealso \code{\link{load_phenoData}}
#' @examples
#' obj = load_meta("~/Desktop/testFile.tsv")
load_meta <-
function(file,sep="\t")
{
	dat2 <- read.table(file,header=FALSE,sep=sep,nrows=1,stringsAsFactors=FALSE);
	subjects <- as.character(dat2[1,-1]);
	classes <-c("character",rep("numeric",length(subjects)));
	dat3 <- read.table(file,header=FALSE,skip=1,sep=sep,colClasses=classes,row.names=1);
	taxa<- rownames(dat3);
	obj <- list(counts=as.data.frame(dat3), taxa=as.data.frame(taxa))
	return(obj);
}
