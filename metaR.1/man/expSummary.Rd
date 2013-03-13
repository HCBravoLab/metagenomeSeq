\name{expSummary}
\Rdversion{1.0}
\docType{methods}
\alias{expSummary,MRexperiment-method}
\alias{expSummary}
\title{
    Access MRexperiment object experiment data
}
\description{
    The expSummary vectors represent the column (sample specific) sums of features, i.e. the total number of reads for a sample, libSize and also the normalization factors, normFactor.
}
\usage{
\S4method{expSummary}{MRexperiment}(obj)
}
\arguments{
  \item{obj}{a \code{MRexperiment} object.}
}

\examples{
data(mouseData)
expSummary(mouseData)
}
\author{
   Joseph N. Paulson, jpaulson@umiacs.umd.edu
}
