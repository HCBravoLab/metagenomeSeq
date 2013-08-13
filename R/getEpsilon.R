#' Calculate the relative difference between iterations of the negative
#' log-likelihoods.
#' 
#' Maximum-likelihood estimates are approximated using the EM algorithm where
#' we treat mixture membership $delta_ij$ = 1 if $y_ij$ is generated from the
#' zero point mass as latent indicator variables. The log-likelihood in this
#' extended model is $(1-delta_ij) log f_count(y;mu_i,sigma_i^2 )+delta_ij log
#' pi_j(s_j)+(1-delta_ij)log (1-pi_j (sj))$. The responsibilities are defined
#' as $z_ij = pr(delta_ij=1 | data)$.
#' 
#' 
#' @param nll Vector of size M with the current negative log-likelihoods.
#' @param nllOld Vector of size M with the previous iterations negative
#' log-likelihoods.
#' @return Vector of size M of the relative differences between the previous
#' and current iteration nll.
#' @seealso \code{\link{fitZig}}
getEpsilon <-
function(nll, nllOld){
	eps=(nllOld-nll)/nllOld
	ifelse(!is.finite(nllOld), Inf, eps)
}
