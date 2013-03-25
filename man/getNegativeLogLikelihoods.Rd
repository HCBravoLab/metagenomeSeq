\name{getNegativeLogLikelihoods}
\alias{getNegativeLogLikelihoods}
\title{Calculate the negative log-likelihoods for the various features given the residuals.}
\usage{
  getNegativeLogLikelihoods(z, countResiduals,
    zeroResiduals)
}
\arguments{
  \item{z}{Matrix (m x n) of estimate responsibilities
  (probabilities that a count comes from a spike
  distribution at 0).}

  \item{countResiduals}{Residuals from the count model.}

  \item{zeroResiduals}{Residuals from the zero model.}
}
\value{
  Vector of size M of the negative log-likelihoods for the
  various features.
}
\description{
  Maximum-likelihood estimates are approximated using the
  EM algorithm where we treat mixture membership
  $delta_{ij}$ = 1 if $y_{ij}$ is generated from the zero
  point mass as latent indicator variables. The
  log-likelihood in this extended model is $(1-delta_{ij})
  log f_{count}(y;mu_i,sigma_i^2 )+delta_{ij} log
  pi_j(s_j)+(1-delta_{ij})log (1-pi_j (sj))$. The
  responsibilities are defined as $z_{ij} =
  pr(delta_{ij}=1 | data and current values)$.
}
\seealso{
  \code{\link{fitZig}}
}

