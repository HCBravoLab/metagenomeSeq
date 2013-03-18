#' Compute the Maximization step calculation for features still active.
#'
#' Maximization step is solved by weighted least squares. The function also computes counts residuals.
#'
#' Maximum-likelihood estimates are approximated using the EM algorithm where we treat mixture membership $\deta_{ij}$ = 1 if $y_{ij}$
#' is generated from the zero point mass as latent indicator variables. The density is defined as $f_zig(y_{ij} = \pi_j(S_j) \cdot f_{0}(y_{ij}) 
#' +(1-\pi_j (S_j))\cdot f_{count}(y_{ij};\mu_i,\sigma_i^2)$.
#' The log-likelihood in this extended model is
#' $(1−\delta_{ij}) \log f_{count}(y;\mu_i,\sigma_i^2 )+\delta_{ij} \log \pi_j(s_j)+(1−\delta_{ij})\log (1−\pi_j (sj))$.
#' The responsibilities are defined as $z_{ij} = pr(\delta_{ij}=1 | data)$.
#'
#' @param z Matrix (m x n) of estimate responsibilities (probabilities that a count comes from a spike distribution at 0).
#' @param y Matrix (m x n) of count observations.
#' @param mmCount Model matrix for the count distribution.
#' @param stillActive Boolean vector of size M, indicating whether a feature converged or not.
#' @param fit2 Previous fit of the count model.
#' @return Update matrix (m x n) of estimate responsibilities (probabilities that a count comes from a spike distribution at 0).
#'
#' @name doCountMStep
#' @seealso \code{\link{fitZig}}
doCountMStep <-
function(z, y, mmCount, stillActive,fit2=NULL){

	if (is.null(fit2)){
	fit=limma::lmFit(y[stillActive,],mmCount,weights = (1-z[stillActive,]))
	countCoef = fit$coefficients
	countMu=tcrossprod(countCoef, mmCount)
	residuals=sweep((y[stillActive,,drop=FALSE]-countMu),1,fit$sigma,"/")

	dat = list(fit = fit, residuals = residuals)
	return(dat)
	} else {

	residuals = fit2$residuals
	fit2 = fit2$fit

	fit=limma::lmFit(y[stillActive,,drop=FALSE],mmCount,weights = (1-z[stillActive,,drop=FALSE]))

	
	fit2$coefficients[stillActive,] = fit$coefficients
	fit2$stdev.unscaled[stillActive,]=fit$stdev.unscaled
	fit2$sigma[stillActive] = fit$sigma
	fit2$Amean[stillActive] = fit$Amean
	fit2$df[stillActive]    = fit$df
	fit2$df.residual[stillActive]    = fit$df.residual

	countCoef = fit$coefficients
	countMu=tcrossprod(countCoef, mmCount)
	r=sweep((y[stillActive,,drop=FALSE]-countMu),1,fit$sigma,"/")
	residuals[stillActive,]=r

	dat = list(fit = fit2, residuals=residuals)

	return(dat)
	}
}
