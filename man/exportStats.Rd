\name{exportStats}
\alias{exportStats}
\title{Various statistics of the count data.}
\usage{
  exportStats(obj, p = cumNormStat(obj),
    output = "~/Desktop/res.stats.tsv")
}
\arguments{
  \item{obj}{A MRexperiment object with count data.}

  \item{p}{Quantile value to calculate the scaling factor
  and quantiles for the various samples.}

  \item{output}{Output file name.}
}
\value{
  None.
}
\description{
  A matrix of values for each sample. The matrix consists
  of sample ids, the sample scaling factor, quantile value,
  the number identified features, and library size (depth
  of coverage).
}
\examples{
# see vignette
}
\seealso{
  \code{\link{cumNorm}} \code{\link{quantile}}
}

