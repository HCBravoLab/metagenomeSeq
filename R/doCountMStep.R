#' Compute the Maximization step calculation for features still active.
#' 
#' Maximization step is solved by weighted least squares.  The function also
#' computes counts residuals.
#' 
#' Maximum-likelihood estimates are approximated using the EM algorithm where
#' we treat mixture membership $delta_ij$ = 1 if $y_ij$ is generated from the
#' zero point mass as latent indicator variables. The density is defined as
#' $f_zig(y_ij = pi_j(S_j)*f_0(y_ij) +(1-pi_j (S_j)) *
#' f_count(y_ij;mu_i,sigma_i^2)$. The log-likelihood in this extended model is
#' $(1-delta_ij) log f_count(y;mu_i,sigma_i^2 )+delta_ij log
#' pi_j(s_j)+(1-delta_ij)log (1-pi_j (s_j))$. The responsibilities are defined
#' as $z_ij = pr(delta_ij=1 | data)$.
#' 
#' @param z Matrix (m x n) of estimate responsibilities (probabilities that a
#' count comes from a spike distribution at 0).
#' @param y Matrix (m x n) of count observations.
#' @param mmCount Model matrix for the count distribution.
#' @param stillActive Boolean vector of size M, indicating whether a feature
#' converged or not.
#' @param fit2 Previous fit of the count model.
#' @param dfMethod Either 'default' or 'modified' (by responsibilities)
#' @return Update matrix (m x n) of estimate responsibilities (probabilities
#' that a count comes from a spike distribution at 0).
#' @seealso \code{\link{fitZig}}
doCountMStep <- function(z, y, mmCount, stillActive,fit2=NULL,dfMethod="modified"){

	if (is.null(fit2)){
	  fit=limma::lmFit(y,mmCount,weights = (1-z))
	  if(dfMethod=="modified"){
	    df = rowSums(1-z) - ncol(mmCount)
	    fit$df = df
	    fit$df.residual = df
	  }
	  countCoef = fit$coefficients
	  countMu=tcrossprod(countCoef, mmCount)
	  residuals=sweep((y-countMu),1,fit$sigma,"/")
	  dat = list(fit = fit, residuals = residuals)
	  return(dat)
		# fit=limma::lmFit(y[stillActive,],mmCount,weights = (1-z[stillActive,]))
		# if(dfMethod=="modified"){
		# 	df = rowSums(1-z[stillActive,,drop=FALSE]) - ncol(mmCount)
		# 	fit$df[stillActive] = df
		# 	fit$df.residual[stillActive] = df
		# }
		# countCoef = fit$coefficients
		# countMu=tcrossprod(countCoef, mmCount)
		# residuals=sweep((y[stillActive,,drop=FALSE]-countMu),1,fit$sigma,"/")
		# dat = list(fit = fit, residuals = residuals)
		# return(dat)
	} else {

		residuals = fit2$residuals
		fit2 = fit2$fit

		fit=limma::lmFit(y[stillActive,,drop=FALSE],mmCount,weights = (1-z[stillActive,,drop=FALSE]))
		
		fit2$coefficients[stillActive,] = fit$coefficients
		fit2$stdev.unscaled[stillActive,]=fit$stdev.unscaled
		fit2$sigma[stillActive] = fit$sigma
		fit2$Amean[stillActive] = fit$Amean

		if(dfMethod=="modified"){
			df = rowSums(1-z[stillActive,,drop=FALSE]) - ncol(mmCount)
			fit$df = df
			fit$df.residual = df
		}
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
