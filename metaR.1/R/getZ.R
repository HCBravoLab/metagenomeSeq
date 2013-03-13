#' Calculate the current Z estimate responsibilities (posterior probabilities)
#'
#' @param z Matrix (m x n) of estimate responsibilities (probabilities that a count comes from a spike distribution at 0).
#' @param zUsed Matrix (m x n) of estimate responsibilities (probabilities that a count comes from a spike distribution at 0) that are actually used (following convergence). 
#' @param stillActive A vector of size M booleans saying if a feature is still active or not.
#' @param nll Vector of size M with the current negative log-likelihoods.
#' @param nllUSED Vector of size M with the converged negative log-likelihoods.
#' @return A list of updated zUsed and nllUSED.
#'
#' @name getZ
#' @seealso \code{\link{fitZig}}
getZ <-
function(z,zUsed,stillActive,nll,nllUSED){

	nllUSED[stillActive] = nll[stillActive]
	k =which(nll< (nllUSED))
	if(length(k)>0){
		zUsed[k,]=z[k,]
		nllUSED[k] = nll[k]
	}
	zUsed[stillActive,] = z[stillActive,]
	dat = list(zUsed = zUsed,nllUSED = nllUSED)
	return(dat);
}
