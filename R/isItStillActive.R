#' Function to determine if a feature is still active.
#'
#' In the Expectation Maximization routine features posterior probabilities routinely converge based on a tolerance threshold. This function checks
#' whether or not the feature's negative log-likelihood (measure of the fit) has changed or not.
#'
#' @param eps Vector of size M (features) representing the relative difference between the new nll and old nll.
#' @param tol The threshold tolerance for the difference
#' @param stillActive A vector of size M booleans saying if a feature is still active or not.
#' @param stillActiveNLL A vector of size M recording the negative log-likelihoods of the various features, updated for those still active.
#' @param nll Vector of size M with the current negative log-likelihoods.
#' @return None.
#'
#' @name isItStillActive
#' @seealso \code{\link{fitZig}}
#'
isItStillActive <-
function(eps, tol,stillActive,stillActiveNLL,nll){
	stillActive[stillActive]=!is.finite(eps[stillActive]) | eps[stillActive]>tol
	stillActive[which(is.na(eps))]=FALSE

	stillActiveNLL[stillActive]=nll[stillActive]
	dat = list(stillActive=stillActive,stillActiveNLL = stillActiveNLL)
	return(dat)
}
