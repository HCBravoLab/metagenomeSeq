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
	dat2 <- read.table(file,header=FALSE,sep="\t");
    # load names 
	subjects <- array(0,dim=c(ncol(dat2)-1));
	for(i in 1:length(subjects)) {
		subjects[i] <- as.character(dat2[1,i+1]);
	}
	
	classes <-c("character",rep("numeric",length(subjects)));
	dat3 <- read.table(file,header=FALSE,skip=1,sep=sep,colClasses=classes);
	
	taxa<- as.matrix(dat3[,1]);

	matrix <- array(0, dim=c(length(taxa),length(subjects)));
	for(i in (1:length(subjects))){
		matrix[,i] <- as.numeric(dat3[,i+1]);
	}	
	
	colnames(matrix) = subjects;
    rownames(matrix) = taxa;
	obj <- list(counts=as.data.frame(matrix), taxa=as.data.frame(taxa))
	
	return(obj);
}
