#' Compute the value of the count density function from the count model residuals.
#'
#' Calculate density values from a normal: $f(x) = 1/(\sqrt (2 \pi ) \sigma ) e^-((x - \mu )^2/(2 \sigma^2))$.
#' Maximum-likelihood estimates are approximated using the EM algorithm where we treat mixture membership $\deta_{ij}$ = 1 if $y_{ij}$
#' is generated from the zero point mass as latent indicator variables. The density is defined as $f_zig(y_{ij} = \pi_j(S_j) \cdot f_{0}(y_{ij}) 
#' +(1-\pi_j (S_j))\cdot f_{count}(y_{ij};\mu_i,\sigma_i^2)$.
#' The log-likelihood in this extended model is
#' $(1−\delta_{ij}) \log f_{count}(y;\mu_i,\sigma_i^2 )+\delta_{ij} \log \pi_j(s_j)+(1−\delta_{ij})\log (1−\pi_j (sj))$.
#' The responsibilities are defined as $z_{ij} = pr(\delta_{ij}=1 | data)$.
#'
#' @param residuals Residuals from the count model.
#' @param log Whether or not we are calculating from a log-normal distribution.
#' @return Density values from the count model residuals.
#'
#' @name getCountDensity
#' @seealso \code{\link{fitZig}}
getCountDensity <-
function(residuals, log=FALSE){
	dnorm(residuals,log=log)
}
