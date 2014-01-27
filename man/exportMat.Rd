\name{exportMat}
\alias{exportMat}
\alias{exportMatrix}
\title{Export the normalized MRexperiment dataset as a matrix.}
\usage{
exportMat(obj, log = TRUE, norm = TRUE, sep = "\\t",
  output = "~/Desktop/matrix.tsv")
}
\arguments{
  \item{obj}{A MRexperiment object or count matrix.}

  \item{log}{Whether or not to log transform the counts -
  if MRexperiment object.}

  \item{norm}{Whether or not to normalize the counts - if
  MRexperiment object.}

  \item{sep}{Separator for writing out the count matrix}

  \item{output}{Output file name}
}
\value{
NA
}
\description{
This function allows the user to take a dataset of counts
and output the dataset to the user's workspace as a
tab-delimited file, etc.
}
\examples{
# see vignette
}
\seealso{
\code{\link{cumNorm}}
}

