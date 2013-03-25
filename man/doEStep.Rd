\name{doEStep}
\alias{doEStep}
\title{Compute the Expectation step.}
\usage{
  doEStep(countResiduals, zeroResiduals, zeroIndices)
}
\arguments{
  \item{countResiduals}{Residuals from the count model.}

  \item{zeroResiduals}{Residuals from the zero model.}

  \item{zeroIndices}{Index (matrix m x n) of counts that
  are zero/non-zero.}
}
\value{
  Updated matrix (m x n) of estimate responsibilities
  (probabilities that a count comes from a spike
  distribution at 0).
}
\description{
  Estimates the responsibilities $z_{ij} = frac{pi_j
  cdot I_{0}(y_{ij}}{pi_j cdot I_{0}(y_{ij} + (1-pi_j)
  cdot f_{count}(y_{ij}}
}
\details{
  Maximum-likelihood estimates are approximated using the
  EM algorithm where we treat mixture membership
  $delta_{ij}$ = 1 if $y_{ij}$ is generated from the zero
  point mass as latent indicator variables. The density is
  defined as $f_zig(y_{ij} = pi_j(S_j) cdot f_{0}(y_{ij})
  +(1-pi_j (S_j))cdot
  f_{count}(y_{ij};mu_i,sigma_i^2)$. The log-likelihood
  in this extended model is $(1-delta_{ij}) log
  f_{count}(y;mu_i,sigma_i^2 )+delta_{ij} log
  pi_j(s_j)+(1-delta_{ij})log (1-pi_j (sj))$. The
  responsibilities are defined as $z_{ij} =
  pr(delta_{ij}=1 | data)$.
}
\seealso{
  \code{\link{fitZig}}
}

