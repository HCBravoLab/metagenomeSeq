\name{normFactors}
\Rdversion{1.0}
\docType{methods}
\alias{normFactors,MRexperiment-method}
\alias{normFactors}
\title{
    Access the normalization factors in a MRexperiment object
}
\description{
  Function to access the scaling factors, aka the normalization factors, of samples in a
    MRexperiment object.
}

\usage{
\S4method{normFactors}{MRexperiment}(obj)
}
\arguments{
  \item{obj}{a \code{MRexperiment} object.}
}
\author{
   Joseph N. Paulson, jpaulson@umiacs.umd.edu
}
\examples{
   data(lungData)
   head(normFactors(lungData))
}   
   
