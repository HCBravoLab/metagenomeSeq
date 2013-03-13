\name{libSize}
\Rdversion{1.0}
\docType{methods}
\alias{libSize,MRexperiment-method}
\alias{libSize}
\title{
    Access sample depth of coverage from MRexperiment object
}
\description{
    The libSize vector  represents the column (sample specific) sums of features, i.e. the total number of reads for a sample. It is used by \code{\link{fitZig}}. 
}
\usage{
\S4method{libSize}{MRexperiment}(obj)
}
\arguments{
  \item{obj}{a \code{MRexperiment} object.}
}

\examples{
   data(lungData)
   head(libSize(lungData))
}
\author{
   Joseph N. Paulson, jpaulson@umiacs.umd.edu
}
