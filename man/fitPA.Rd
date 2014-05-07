\name{fitPA}
\alias{fitPA}
\title{Wrapper to run fisher's test on presence/absence of a feature.}
\usage{
fitPA(obj, cl, thres = 0, parallel = FALSE, cores = 2, ...)
}
\arguments{
  \item{obj}{A MRexperiment object with a count matrix, or
  a simple count matrix.}

  \item{cl}{Group comparison}

  \item{thres}{Threshold for defining presence/absence.}

  \item{parallel}{Use multiple cores?}

  \item{cores}{Number of cores to use.}

  \item{...}{Extra options for makeCluster}
}
\value{
NA
}
\description{
This function returns a data frame of p-values, odds
ratios, lower and upper confidence limits for every row of
a matrix.
}
\examples{
data(lungData)
k = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-k]
lungTrim = lungTrim[-which(rowSums(MRcounts(lungTrim)>0)<20),]
res = fitPA(lungTrim,pData(lungTrim)$SmokingStatus);
head(res)
}
\seealso{
\code{\link{cumNorm}} \code{\link{fitZig}}
}

