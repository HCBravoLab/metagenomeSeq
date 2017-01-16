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
#' @param dampening_limit <numeric> z values used will be within param and 1-param (default: 1e-8)
#' @param per_feature <logical> fit zero model per feature independently (default: FALSE)
#' 
#' @return List of the zero fit (zero mean model) coefficients, variance -
#' scale parameter (scalar), and normalized residuals of length
#' sum(zeroIndices).
#' @seealso \code{\link{fitZig}}
doZeroMStep <- function(z, zeroIndices, mmZero, dampening_limit=1e-8, per_feature=FALSE) {
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
    .dampen(tmp, dampening_limit)
  })
  
  zeroLM <- lm.fit(mmZero, qlogis(pi))
  zeroCoef <- zeroLM$coef
  
  r <- zeroLM$residuals
  sigma <- sd(r) + (1e-3)
  
  list(zeroLM=zeroLM, zeroCoef=zeroCoef, sigma=sigma, residuals=r/sigma)
}

.dampen <- function(y, dampening_limit) {
  pmax(dampening_limit, pmin(1-dampening_limit, y))
}

.doZeroMStep_perFeature <- function(z, zeroIndices, mmZero, dampening_limit=1e-8) {
  nfeatures <- nrow(z)
  zeroLMs <- vector("list", nfeatures)
  zeroCoef <- matrix(NA, nrow=nfeatures, ncol=ncol(mmZero))
  r <- matrix(NA, nrow=nfeatures, ncol=ncol(z))
  sigma <- vector("numeric", nfeatures)

  for (i in seq_len(nfeatures)) {
    res <- .doZeroMStep_oneFeature(z[i, ], zeroIndices[i, ], mmZero, dampening_limit)
    zeroLMs[[i]] <- res$fit
    zeroCoef[i, ] <- res$coef
    r[i, ] <- res$residuals
    sigma[i] <- res$sigma
  }
  
  list(zeroLM=zeroLMs, zeroCoef=zeroCoef, sigma=sigma, residuals=r)
}

.doZeroMStep_oneFeature <- function(z, zeroIndices, mmZero, dampening_limit) {
  y <- z[zeroIndices]
  y <- .dampen(y, dampening_limit)
  y <- qlogis(y)
  mmZero <- mmZero[zeroIndices,,drop=FALSE]
  fit <- lm.fit(mmZero, y)
  coef <- fit$coef
  r <- vector("numeric", length(z))
  r[zeroIndices] <- fit$residuals
  r[!zeroIndices] <- NA
  sigma <- sd(r, na.rm=TRUE)
  if(sigma == 0) sigma = 1
  r <- r / sigma
  list(fit=fit, coef=coef, residuals=r, sigma=sigma)
}
