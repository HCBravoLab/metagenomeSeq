\docType{methods}
\name{libSize}
\alias{libSize}
\alias{libSize,MRexperiment-method}
\title{Access sample depth of coverage from MRexperiment object}
\usage{
libSize(obj)
}
\arguments{
  \item{obj}{a \code{MRexperiment} object.}
}
\description{
The libSize vector represents the column (sample specific)
sums of features, i.e. the total number of reads for a
sample or depth of coverage. It is used by
\code{\link{fitZig}}.
}
\examples{
data(lungData)
head(libSize(lungData))
}
\author{
Joseph N. Paulson, jpaulson@umiacs.umd.edu
}

