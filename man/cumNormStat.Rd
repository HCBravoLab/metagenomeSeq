\name{cumNormStat}
\alias{cumNormStat}
\title{Cumulative sum scaling percentile selection}
\usage{
  cumNormStat(obj, qFlag = TRUE, pFlag = FALSE, rel = 0.1,
    ...)
}
\arguments{
  \item{obj}{A list with count data}

  \item{qFlag}{Flag to either calculate the proper
  percentile using a step-wise or triangular approximation
  of the sample count distribution (default step-wise).}

  \item{pFlag}{Plot the median difference quantiles}

  \item{rel}{Cutoff for the relative difference from one
  median difference from the reference to the next}

  \item{...}{Applicable if pFlag == TRUE. Extra plotting
  parameters.}
}
\value{
  Percentile for which to scale data
}
\description{
  Calculates the percentile for which to sum counts up to
  and scale by.
}
\examples{
data(mouseData)
p = round(cumNormStat(mouseData,pFlag=FALSE),digits=2)
}
\seealso{
  \code{\link{fitZig}} \code{\link{cumNorm}}
}

