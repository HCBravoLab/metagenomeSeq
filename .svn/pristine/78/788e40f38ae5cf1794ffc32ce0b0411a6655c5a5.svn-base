\name{getEpsilon}
\alias{getEpsilon}
\title{Calculate the relative difference between iterations of the negative
log-likelihoods.}
\usage{
getEpsilon(nll, nllOld)
}
\arguments{
  \item{nll}{Vector of size M with the current negative
  log-likelihoods.}

  \item{nllOld}{Vector of size M with the previous
  iterations negative log-likelihoods.}
}
\value{
Vector of size M of the relative differences between the
previous and current iteration nll.
}
\description{
Maximum-likelihood estimates are approximated using the EM
algorithm where we treat mixture membership $delta_ij$ = 1
if $y_ij$ is generated from the zero point mass as latent
indicator variables. The log-likelihood in this extended
model is $(1-delta_ij) log f_count(y;mu_i,sigma_i^2
)+delta_ij log pi_j(s_j)+(1-delta_ij)log (1-pi_j (sj))$.
The responsibilities are defined as $z_ij = pr(delta_ij=1 |
data)$.
}
\seealso{
\code{\link{fitZig}}
}

