#' Compute the Expectation step.
#' 
#' Estimates the responsibilities $z_ij = fracpi_j cdot I_0(y_ijpi_j cdot
#' I_0(y_ij + (1-pi_j) cdot f_count(y_ij
#' 
#' Maximum-likelihood estimates are approximated using the EM algorithm where
#' we treat mixture membership $delta_ij$ = 1 if $y_ij$ is generated from the
#' zero point mass as latent indicator variables. The density is defined as
#' $f_zig(y_ij = pi_j(S_j) cdot f_0(y_ij) +(1-pi_j (S_j))cdot
#' f_count(y_ij;mu_i,sigma_i^2)$. The log-likelihood in this extended model is
#' $(1-delta_ij) log f_count(y;mu_i,sigma_i^2 )+delta_ij log
#' pi_j(s_j)+(1-delta_ij)log (1-pi_j (sj))$. The responsibilities are defined
#' as $z_ij = pr(delta_ij=1 | data)$.
#' 
#' @param countResiduals Residuals from the count model.
#' @param zeroResiduals Residuals from the zero model.
#' @param zeroIndices Index (matrix m x n) of counts that are zero/non-zero.
#' @return Updated matrix (m x n) of estimate responsibilities (probabilities
#' that a count comes from a spike distribution at 0).
#' @seealso \code{\link{fitZig}}
doEStep <- function(countResiduals,  zeroResiduals, zeroIndices, per_feature=FALSE) {
  
	pi_prop <- getPi(zeroResiduals)
	
	if (!isTRUE(per_feature)) {
  	w1 <- sweep(zeroIndices, 2, pi_prop, FUN="*")
	} else {
	  w1 <- pi_prop
	}

	countDensity <- getCountDensity(countResiduals)
	
	if (!isTRUE(per_feature)) {
	  w2 <- sweep(countDensity, 2, 1-pi_prop, FUN="*")
	} else {
	  w2 <- countDensity * 1-pi_prop
	}
	
	z <- w1/(w1+w2)
	z[z>1-1e-6] <- 1-1e-6
	z[!zeroIndices] <- 0
	z
}
