#' Load a count dataset associated with a study set up in a Qiime format.
#'
#' Load a matrix of OTUs in Qiime's format
#'
#' @param file Path and filename of the actual data file.
#' @return An object of count data.
#'
#' @name load_metaQ
#' @aliases qiimeLoader
#' @seealso \code{\link{load_meta}} \code{\link{load_phenoData}}
#' @examples
#' obj = load_metaQ("~/Desktop/testFile.tsv")
load_metaQ <- function(file) {	
	dat2 <- read.delim(file,header=FALSE,stringsAsFactors=FALSE,nrows=1,skip=1);
	len = ncol(dat2)
	subjects = as.character(dat2[1,-c(1,len)]);
	classes <-c("character",rep("numeric",(len-2)),"character");
	dat3 <- read.delim(file,header=TRUE,colClasses=classes,skip=1);
	taxa<- dat3[,len];
	taxa<-as.matrix(taxa);
	matrix <- dat3[,-c(1,len)]
	colnames(matrix) = subjects;
	otus = dat3[,1];
	rownames(matrix) = otus;
	obj <- list(counts=as.data.frame(matrix), taxa=as.data.frame(taxa),otus = as.data.frame(otus))
	return(obj);
}
