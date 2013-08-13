\name{MRfisher}
\alias{MRfisher}
\title{Wrapper to run fisher's test on presence/absence of a feature.}
\usage{
  MRfisher(obj, cl, mat = FALSE)
}
\arguments{
  \item{obj}{A MRexperiment object with a count matrix, or
  a simple count matrix.}

  \item{cl}{Group comparison}

  \item{mat}{logical indicating whether obj is a
  MRexperiment object or matrix. Default is a MRexperiment
  object.}
}
\value{
  NA
}
\description{
  This function returns a data frame of p-values, odds
  ratios, lower and upper confidence limits for every row
  of a matrix.
}
\examples{
data(lungData)
k = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-k]
lungTrim = lungTrim[-which(rowSums(MRcounts(lungTrim)>0)<20),]
res = MRfisher(lungTrim,pData(lungTrim)$SmokingStatus);
head(res)
}
\seealso{
  \code{\link{cumNorm}} \code{\link{fitZig}}
}

