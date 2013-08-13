#' export the normalized eSet dataset as a matrix.
#' 
#' This function allows the user to take a dataset of counts and output the
#' dataset to the user's workspace as a tab-delimited file, etc.
#' 
#' 
#' @aliases exportMatrix exportMat
#' @param mat A matrix of values (normalized, or otherwise)
#' @param output Output file name
#' @return NA
#' @seealso \code{\link{cumNorm}}
#' @examples
#' 
#' # see vignette
#' 
exportMat <-
function(mat,output="~/Desktop/matrix.tsv"){
	matrix = mat;
	
	mat = array(NA,dim=c((nrow(matrix)+1),(ncol(matrix)+1)));
	mat[1,2:ncol(mat)] = colnames(matrix);
	mat[2:nrow(mat),2:ncol(mat)] = matrix;
    mat[2:nrow(mat),1] = rownames(matrix);
    mat[1,1] = "Taxa and Samples"
	write(t(mat),file=output,sep="\t",ncolumns=ncol(mat))	
}
