\name{cumNormMat}
\alias{cumNormMat}
\title{Cumulative sum scaling factors.}
\usage{
  cumNormMat(obj, p = cumNormStat(obj))
}
\arguments{
  \item{obj}{A MRexperiment object.}

  \item{p}{The pth quantile.}
}
\value{
  Returns a matrix normalized by scaling counts up to and including the pth quantile.
}
\description{
  Calculates each column's quantile and calculates the sum
  up to and including that quantile.
}
\examples{
data(mouseData)
head(cumNormMat(mouseData))
}
\seealso{
  \code{\link{fitZig}} \code{\link{cumNorm}}
}

