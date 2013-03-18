\name{newMRexperiment}
\Rdversion{1.0}
\alias{newMRexperiment}
\title{
Create a MRexperiment object
}
\description{
This function creates a MRexperiment object from a matrix or data frame of count data.
}
\usage{
newMRexperiment(counts, phenoData=NULL, featureData=NULL,libSize=NULL, normFactors=NULL)
}
\arguments{
  \item{counts}{
     A matrix or data frame of count data. The count data is representative of the number of reads annotated for a feature (be it gene, OTU, species, etc). Rows should correspond to features and columns to samples.
}
  \item{phenoData}{
      An AnnotatedDataFrame with pertinent sample information.
}
  \item{featureData}{
      An AnnotatedDataFrame with pertinent feature information.
}
  \item{libSize}{
      libSize, library size, is the total number of reads for a particular sample.
}
  \item{normFactors}{
      normFactors, the normalization factors used in either the model or as scaling factors of sample counts for each particular sample.
}
}
\details{
   See  \code{\link{MRexperiment-class}} and \code{eSet} (from the Biobase package) for the meaning of the various slots.
}
\value{
   an object of class MRexperiment
}
\author{
   Joseph N Paulson, jpaulson@umiacs.umd.edu
}
\examples{
cnts = matrix(abs(rnorm(1000)),nc=10)
obj <- newMRexperiment(cnts)
}
