#' export the normalized eSet dataset as a matrix.
#'
#' This function allows the user to take the normalized dataset or counts and output
#' the dataset to the user's workspace as a tab-delimited file, etc.
#'
#' @param jobj An eSet object with count data.
#' @param output Output file name.
#' @return NA
#'
#' @name export_mat
#' @aliases exportMatrix
#' @seealso \code{\link{cumNorm}}
#' @examples
#' export_mat(jobj,output="~/Desktop/normMatrix.tsv");

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
