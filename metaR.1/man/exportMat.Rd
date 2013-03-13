\name{exportMat}
\alias{exportMatrix}
\alias{exportMat}
\title{export the normalized eSet dataset as a matrix.}
\usage{
  exportMat(mat, output = "~/Desktop/matrix.tsv")
}
\arguments{
  \item{mat}{A matrix of values (normalized, or otherwise)}
  \item{output}{Output file name}
}
\value{
  NA
}
\description{
  This function allows the user to take a
  dataset of counts and output the dataset to the user's
  workspace as a tab-delimited file, etc.
}
\examples{
# see vignette
}
\seealso{
  \code{\link{cumNorm}}
}

