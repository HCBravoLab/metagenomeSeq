\name{cumNorm}
\alias{cumNorm}
\title{Cumulative sum scaling factors.}
\usage{
  cumNorm(obj, p = cumNormStat(obj))
}
\arguments{
  \item{obj}{An MRexperiment object.}

  \item{p}{The pth quantile.}
}
\value{
  Vector of the sum up to and including a sample's pth
  quantile
}
\description{
  Calculates each column's quantile and calculates the sum
  up to and including that quantile.
}
\examples{
data(mouseData)
cumNorm(mouseData)
head(normFactors(mouseData))
}
\seealso{
  \code{\link{fitZig}} \code{\link{cumNormStat}}
}

