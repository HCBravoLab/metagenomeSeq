#' Settings for the fitZig function
#'
#' @param tol The tolerance for the difference in negative log likelihood estimates for a feature to remain active.
#' @param maxit The maximum number of iterations for the expectation-maximization algorithm.
#' @param verbose Whether to display iterative step summary statistics or not.
#' @param dfMethod Method to calculate degrees of freedom, either 'default', (non-zero) or 'modified' (weighted by responsibilities).
#' @param pvalMethod Method to calculate p-values, either 'default' (parametric) or 'bootstrap'.
#' @param per_feature_zeroModel <logical> Fit a logistic model per feature indepdently (default: FALSE)
#' @return The value for the tolerance, maximum no. of iterations, and the verbose warning.
#' @note \code{\link{fitZig}} makes use of zigControl.
#'
#' @name zigControl
#' @aliases settings2
#' @seealso \code{\link{fitZig}} \code{\link{cumNorm}} \code{\link{plotOTU}}
#' @examples
#' control =  zigControl(tol=1e-10,maxit=10,verbose=FALSE)
#'
zigControl <- function(tol=1e-4,
                      maxit=10,
                      verbose=TRUE,
                      dfMethod=c("modified", "default"),
                      pvalMethod=c("default", "bootstrap"),
                      per_feature_zeroModel=FALSE) {

	dfMethod <- match.arg(dfMethod)
	pvalMethod <- match.arg(pvalMethod)
	
	set <-list(tol=tol,
	           maxit=maxit,
	           verbose=verbose,
	           dfMethod=dfMethod,
	           pvalMethod=pvalMethod,
	           per_feature_zeroModel=per_feature_zeroModel)
	return(set)	
}
