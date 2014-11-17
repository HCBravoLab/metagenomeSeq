#' Settings for the fitZig function
#'
#' @param tol The tolerance for the difference in negative log likelihood estimates for a feature to remain active.
#' @param maxit The maximum number of iterations for the expectation-maximization algorithm.
#' @param verbose Whether to display iterative step summary statistics or not.
#' @param dfmethod Either 'default' or 'modified' (by responsibilities)
#' @return The value for the tolerance, maximum no. of iterations, and the verbose warning.
#' @note \code{\link{fitZig}} makes use of zigControl.
#'
#' @name zigControl
#' @aliases settings2
#' @seealso \code{\link{fitZig}} \code{\link{cumNorm}} \code{\link{plotOTU}}
#' @examples
#' control =  zigControl(tol=1e-10,maxit=10,verbose=FALSE)
#'
zigControl <-function(tol=1e-4,maxit=10,verbose=TRUE,dfmethod="default"){
	set <-list(tol=tol,maxit=maxit,verbose=verbose,dfmethod=dfmethod);
	return(set)	
}
