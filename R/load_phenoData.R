#' Load a clinical/phenotypic dataset associated with a study.
#' 
#' Load a matrix of metadata associated with a study.
#' 
#' 
#' @aliases load_phenoData phenoData
#' @param file Path and filename of the actual clinical file.
#' @param tran Boolean. If the covariates are along the columns and samples
#' along the rows, then tran should equal TRUE.
#' @param sep The separator for the file.
#' @return The metadata as a dataframe.
#' @seealso \code{\link{load_meta}}
#' @examples
#' 
#' # see vignette
#' 
load_phenoData <-
function(file,tran=FALSE,sep="\t")
{
	dat2 <- read.table(file,header=FALSE,sep=sep);

# no. of subjects 
	subjects <- array(0,dim=c(ncol(dat2)-1));
	for(i in 1:length(subjects)) {
		subjects[i] <- as.character(dat2[1,i+1]);
	}
# no. of rows
	rows <- nrow(dat2);
	
# load remaining counts
	matrix <- array(NA, dim=c(length(subjects),rows-1));
	covar = array(NA,dim=c(rows-1,1)); 
	
	for(i in 1:(rows)-1){
		for(j in 1:(length(subjects))){ 
			matrix[j,i] <- as.character(dat2[i+1,j+1]);
		}
		covar[i] = as.character(dat2[i+1,1]);
	}  
	
		
	phenoData<-as.data.frame(matrix);
	
	colnames(phenoData) = covar;
	if(length(unique(subjects))==(length(subjects))){
		rownames(phenoData) = subjects;
	}
    if(tran==TRUE){
        phenoData = as.data.frame(t(phenoData))
	}
	return(phenoData);
}
