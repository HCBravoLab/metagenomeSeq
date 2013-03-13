#' Calculate the negative log-likelihoods for the various features given the residuals.
#'
#' Maximum-likelihood estimates are approximated using the EM algorithm where we treat mixture membership $\deta_{ij}$ = 1 if $y_{ij}$
#' is generated from the zero point mass as latent indicator variables. The log-likelihood in this extended model is
#' $(1−\delta_{ij}) \log f_{count}(y;\mu_i,\sigma_i^2 )+\delta_{ij} \log \pi_j(s_j)+(1−\delta_{ij})\log (1−\pi_j (sj))$.
#' The responsibilities are defined as $z_{ij} = pr(\delta_{ij}=1 | data and current  values)$.
#'
#' @param z Matrix (m x n) of estimate responsibilities (probabilities that a count comes from a spike distribution at 0).
#' @param countResiduals Residuals from the count model.
#' @param zeroResiduals Residuals from the zero model.
#' @return Vector of size M of the negative log-likelihoods for the various features.
#'
#' @name getNegativeLogLikelihoods
#' @seealso \code{\link{fitZig}}
getNegativeLogLikelihoods <-
function(z, countResiduals, zeroResiduals){
	pi=getPi(zeroResiduals)
	countDensity=getCountDensity(countResiduals, log=TRUE)
	res=(1-z) * countDensity
	res=res+sweep(z, 2, log(pi), FUN="*")
	res=res+sweep(1-z,2,log(1-pi), FUN="*")
	-rowSums(res)
}
