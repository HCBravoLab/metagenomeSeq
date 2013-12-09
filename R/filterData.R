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
    if(class(obj)=="MRexperiment"){
        mat = MRcounts(obj,norm=FALSE,log=FALSE)>0
    } else if(class(obj) == "matrix") {
        mat = obj>0
    } else {
        stop("Object needs to be either a MRexperiment object or matrix.")
    }
	cols = which(colSums(obj)>=depth)
	rows = which(rowSums(mat[,cols])>=present)
	return(obj[rows,cols])
}
