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
load_metaQ <-
function(file)
{	
	dat2 <- read.table(file,header=FALSE,sep="\t");
# load names 
	subjects <- array(0,dim=c(ncol(dat2)-1));
	for(i in 1:length(subjects)) {
		subjects[i] <- as.character(dat2[1,i+1]);
	}
    classes <-c("character",rep("numeric",(length(subjects)-1)),"character");
	dat3 <- read.table(file,header=TRUE,sep="\t",colClasses=classes);
	
	taxa<- dat3[,1+length(subjects)];
	taxa<-as.matrix(taxa);
	
	matrix <- array(0, dim=c(length(taxa),(length(subjects)-1)));
	for(i in (1:(length(subjects)-1))){
		matrix[,i] <- as.numeric(dat3[,i+1]);
	}	
	
	colnames(matrix) = subjects[-length(subjects)];
	obj <- list(counts=matrix, taxa=as.data.frame(taxa),otus = as.data.frame(dat3[,1]))
	
	return(obj);
}
