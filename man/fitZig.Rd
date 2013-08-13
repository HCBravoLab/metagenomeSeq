\name{fitZig}
\alias{fitZig}
\title{Computes the weighted fold-change estimates and t-statistics.}
\usage{
  fitZig(obj, mod, zeroMod = NULL, useS95offset = TRUE,
    control = zigControl())
}
\arguments{
  \item{obj}{A MRexperiment object with count data.}

  \item{mod}{The model for the count distribution.}

  \item{zeroMod}{The zero model, the model to account for
  the change in the number of OTUs observed as a linear
  effect of the depth of coverage.}

  \item{useS95offset}{Boolean, whether to include the
  default scaling parameters in the model or not.}

  \item{control}{The settings for fitZig.}
}
\value{
  The fits, posterior probabilities, posterior
  probabilities used at time of convergence for each
  feature, ebayes (limma object) fit, among other data.
}
\description{
  Wrapper to actually run the Expectation-maximization
  algorithm and estimate $f_count$ fits.
  Maximum-likelihood estimates are approximated using the
  EM algorithm where we treat mixture membership $delta_ij
  = 1$ if $y_ij$ is generated from the zero point mass as
  latent indicator variables. The density is defined as
  $f_zig(y_ij = pi_j(S_j)*f_0(y_ij) +(1-pi_j (S_j)) *
  f_count(y_ij; mu_i, sigma_i^2)$. The log-likelihood in
  this extended model is: $(1-delta_ij) log
  f_count(y;mu_i,sigma_i^2 )+delta_ij log
  pi_j(s_j)+(1-delta_ij) log (1-pi_j (s_j))$. The
  responsibilities are defined as $z_ij = pr(delta_ij=1 |
  data)$.
}
\examples{
data(lungData)
k = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-k]
k = which(rowSums(MRcounts(lungTrim)>0)<30)
cumNorm(lungTrim)
lungTrim = lungTrim[-k,]
smokingStatus = pData(lungTrim)$SmokingStatus
mod = model.matrix(~smokingStatus)
settings = zigControl(maxit=1,verbose=FALSE)
fit = fitZig(obj = lungTrim,mod=mod,control=settings)
}
\seealso{
  \code{\link{cumNorm}} \code{\link{zigControl}}
}

