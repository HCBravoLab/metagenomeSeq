#' Compute the Expectation step.
#'
#' Estimates the responsibilities $z_{ij} = \frac{\pi_j \cdot I_{0}(y_{ij}}{\pi_j \cdot I_{0}(y_{ij} + (1-\pi_j) \cdot f_{count}(y_{ij}}
#'
#' Maximum-likelihood estimates are approximated using the EM algorithm where we treat mixture membership $\deta_{ij}$ = 1 if $y_{ij}$
#' is generated from the zero point mass as latent indicator variables. The density is defined as $f_zig(y_{ij} = \pi_j(S_j) \cdot f_{0}(y_{ij}) 
#' +(1-\pi_j (S_j))\cdot f_{count}(y_{ij};\mu_i,\sigma_i^2)$.
#' The log-likelihood in this extended model is
#' $(1−\delta_{ij}) \log f_{count}(y;\mu_i,\sigma_i^2 )+\delta_{ij} \log \pi_j(s_j)+(1−\delta_{ij})\log (1−\pi_j (sj))$.
#' The responsibilities are defined as $z_{ij} = pr(\delta_{ij}=1 | data)$.
#'
#' @param countResiduals Residuals from the count model.
#' @param zeroResiduals Residuals from the zero model.
#' @param zeroIndices Index (matrix m x n) of counts that are zero/non-zero.
#' @return Updated matrix (m x n) of estimate responsibilities (probabilities that a count comes from a spike distribution at 0).
#'
#' @name doEStep
#' @seealso \code{\link{fitZig}}
doEStep <-
function(countResiduals,  zeroResiduals, zeroIndices)
{
	pi_prop=getPi(zeroResiduals)
	w1=sweep(zeroIndices, 2, pi_prop, FUN="*")

	countDensity=getCountDensity(countResiduals)
	w2=sweep(countDensity, 2, 1-pi_prop, FUN="*")
	z=w1/(w1+w2)
	z[!zeroIndices]=0
	z
}
