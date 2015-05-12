#' Filter datasets according to no. features present in features with at least a certain depth.
#' 
#' Filter the data based on the number of present features after filtering samples by depth of coverage.
#' There are many ways to filter the object, this is just one way.
#' 
#' @param obj A MRexperiment object or count matrix.
#' @param present Features with at least 'present' postive samples.
#' @param depth Sampls with at least this much depth of coverage
#' @return A MRexperiment object.
#' @export
#' @examples
#' 
#' data(mouseData)
#' filterData(mouseData)
#' 
filterData <- function(obj,present=1,depth=1000){
    mat = returnAppropriateObj(obj,norm=FALSE,log=FALSE)>0
	cols = which(colSums(MRcounts(obj))>=depth)
	rows = which(rowSums(mat[,cols])>=present)
	return(obj[rows,cols])
}
