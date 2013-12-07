\name{isItStillActive}
\alias{isItStillActive}
\title{Function to determine if a feature is still active.}
\usage{
isItStillActive(eps, tol, stillActive, stillActiveNLL, nll)
}
\arguments{
  \item{eps}{Vector of size M (features) representing the
  relative difference between the new nll and old nll.}

  \item{tol}{The threshold tolerance for the difference}

  \item{stillActive}{A vector of size M booleans saying if
  a feature is still active or not.}

  \item{stillActiveNLL}{A vector of size M recording the
  negative log-likelihoods of the various features, updated
  for those still active.}

  \item{nll}{Vector of size M with the current negative
  log-likelihoods.}
}
\value{
None.
}
\description{
In the Expectation Maximization routine features posterior
probabilities routinely converge based on a tolerance
threshold. This function checks whether or not the
feature's negative log-likelihood (measure of the fit) has
changed or not.
}
\seealso{
\code{\link{fitZig}}
}

