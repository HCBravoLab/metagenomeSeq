\name{cumNormStatFast}
\alias{cumNormStatFast}
\title{Cumulative sum scaling percentile selection}
\usage{
  cumNormStatFast(obj, pFlag = FALSE, rel = 0.1, ...)
}
\arguments{
  \item{obj}{A list with count data}

  \item{pFlag}{Plot the median difference quantiles}

  \item{rel}{Cutoff for the relative difference from one
  median difference from the reference to the next}

  \item{...}{Applicable if pFlag == TRUE. Additional
  plotting parameters.}
}
\value{
  Percentile for which to scale data
}
\description{
  Calculates the percentile for which to sum counts up to
  and scale by. Faster version than available in
  cumNormStat.
}
\examples{
data(mouseData)
p = round(cumNormStatFast(mouseData,pFlag=FALSE),digits=2)
}
\seealso{
  \code{\link{fitZig}} \code{\link{cumNorm}}
  \code{\link{cumNormStat}}
}

