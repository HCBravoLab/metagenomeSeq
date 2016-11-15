#' Compute the zero Maximization step.
#' 
#' Performs Maximization step calculation for the mixture components. Uses
#' least squares to fit the parameters of the mean of the logistic
#' distribution. $$ pi_j = sum_i^M frac1Mz_ij $$ Maximum-likelihood estimates
#' are approximated using the EM algorithm where we treat mixture membership
#' $delta_ij$ = 1 if $y_ij$ is generated from the zero point mass as latent
#' indicator variables. The density is defined as $f_zig(y_ij = pi_j(S_j) cdot
#' f_0(y_ij) +(1-pi_j (S_j))cdot f_count(y_ij;mu_i,sigma_i^2)$. The
#' log-likelihood in this extended model is $(1-delta_ij) log
#' f_count(y;mu_i,sigma_i^2 )+delta_ij log pi_j(s_j)+(1-delta_ij)log (1-pi_j
#' (sj))$. The responsibilities are defined as $z_ij = pr(delta_ij=1 | data)$.
#' 
#' 
#' @param z Matrix (m x n) of estimate responsibilities (probabilities that a
#' count comes from a spike distribution at 0).
#' @param zeroIndices Index (matrix m x n) of counts that are zero/non-zero.
#' @param mmZero The zero model, the model matrix to account for the change in
#' the number of OTUs observed as a linear effect of the depth of coverage.
#' @return List of the zero fit (zero mean model) coefficients, variance -
#' scale parameter (scalar), and normalized residuals of length
#' sum(zeroIndices).
#' @seealso \code{\link{fitZig}}
doZeroMStep <-
function(z, zeroIndices, mmZero, dampening_limit=1e-8, per_feature=FALSE)
{
  if (isTRUE(per_feature)) {
    .doZeroMStep_perFeature(z, zeroIndices, mmZero, dampening_limit)
  } else {
    .doZeroMStep_combined(z, zeroIndices, mmZero, dampening_limit)
  }
}

.doZeroMStep_combined <- function(z, zeroIndices, mmZero, dampening_limit=1e-8) {
  pi <- sapply(1:ncol(zeroIndices), function(j) {
    if (sum(zeroIndices[,j]) == 0) {
      return(dampening_limit)
    }
    
    tmp <- mean(z[zeroIndices[,j],j], na.rm=TRUE)
    ifelse(tmp <= 1e-8, 
           dampening_limit, 
           ifelse( tmp >= 1-dampening_limit, 
                   1-dampening_limit, 
                   tmp)) 
  })
  
  zeroLM <- lm.fit(mmZero, qlogis(pi))
  zeroCoef <- zeroLM$coef
  
  r <- zeroLM$residuals
  sigma <- sd(r) + (1e-3)
  
  list(zeroLM=zeroLM, zeroCoef=zeroCoef, sigma=sigma, residuals=r/sigma)
}

.doZeroMStep_perFeature <- function(z, zeroIndices, mmZero, dampening_limit=1e-8) {
  stop("not implemented")
}
