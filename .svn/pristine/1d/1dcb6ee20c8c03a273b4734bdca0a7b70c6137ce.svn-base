\name{doCountMStep}
\alias{doCountMStep}
\title{Compute the Maximization step calculation for features still active.}
\usage{
doCountMStep(z, y, mmCount, stillActive, fit2 = NULL)
}
\arguments{
  \item{z}{Matrix (m x n) of estimate responsibilities
  (probabilities that a count comes from a spike
  distribution at 0).}

  \item{y}{Matrix (m x n) of count observations.}

  \item{mmCount}{Model matrix for the count distribution.}

  \item{stillActive}{Boolean vector of size M, indicating
  whether a feature converged or not.}

  \item{fit2}{Previous fit of the count model.}
}
\value{
Update matrix (m x n) of estimate responsibilities
(probabilities that a count comes from a spike distribution
at 0).
}
\description{
Maximization step is solved by weighted least squares.  The
function also computes counts residuals.
}
\details{
Maximum-likelihood estimates are approximated using the EM
algorithm where we treat mixture membership $delta_ij$ = 1
if $y_ij$ is generated from the zero point mass as latent
indicator variables. The density is defined as $f_zig(y_ij
= pi_j(S_j)*f_0(y_ij) +(1-pi_j (S_j)) *
f_count(y_ij;mu_i,sigma_i^2)$. The log-likelihood in this
extended model is $(1-delta_ij) log
f_count(y;mu_i,sigma_i^2 )+delta_ij log
pi_j(s_j)+(1-delta_ij)log (1-pi_j (s_j))$. The
responsibilities are defined as $z_ij = pr(delta_ij=1 |
data)$.
}
\seealso{
\code{\link{fitZig}}
}

