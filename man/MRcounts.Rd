\name{MRcounts}
\Rdversion{1.0}
\docType{methods}
\alias{MRcounts,MRexperiment-method}
\alias{MRcounts}
\title{
   Accessor for the counts slot of a MRexperiment object
}

\description{
   The counts slot holds the raw count data representing (along the rows) the number of reads annotated for a particular feature and (along the columns) the sample.
}

\usage{
\S4method{MRcounts}{MRexperiment}(cnts, norm=FALSE)
}

\arguments{
  \item{cnts}{a \code{MRexperiment} object.}
  \item{norm}{logical indicating whether or not to return normalized counts.}
}

\author{
   Joseph N. Paulson, jpaulson@umiacs.umd.edu
}
\examples{
   data(lungData)
   head(MRcounts(lungData))
}   
   
