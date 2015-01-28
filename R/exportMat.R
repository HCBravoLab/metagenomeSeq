#' Export the normalized MRexperiment dataset as a matrix.
#' 
#' This function allows the user to take a dataset of counts and output the
#' dataset to the user's workspace as a tab-delimited file, etc.
#' 
#' 
#' @aliases exportMatrix exportMat
#' @param obj A MRexperiment object or count matrix.
#' @param log Whether or not to log transform the counts - if MRexperiment object.
#' @param norm Whether or not to normalize the counts - if MRexperiment object.
#' @param sep Separator for writing out the count matrix.
#' @param file Output file name.
#' @return NA
#' @seealso \code{\link{cumNorm}}
#' @examples
#' 
#' # see vignette
#' 
exportMat <-function(obj,log=TRUE,norm=TRUE,sep="\t",file="~/Desktop/matrix.tsv"){
    mat = returnAppropriateObj(obj,norm,log)
	oMat = array(NA,dim=c((nrow(mat)+1),(ncol(mat)+1)));
	oMat[1,2:ncol(oMat)] = colnames(mat);
	oMat[2:nrow(oMat),2:ncol(oMat)] = mat;
    oMat[2:nrow(oMat),1] = rownames(mat);
    oMat[1,1] = "Taxa and Samples";
	write(t(oMat),file=file,sep=sep,ncolumns=ncol(oMat));
}
