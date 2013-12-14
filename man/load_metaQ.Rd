\name{load_metaQ}
\alias{load_metaQ}
\alias{qiimeLoader}
\title{Load a count dataset associated with a study set up in a Qiime format.}
\usage{
load_metaQ(file)
}
\arguments{
  \item{file}{Path and filename of the actual data file.}
}
\value{
An list with 'counts' containing the count data, 'taxa'
containing the otu annotation, and 'otus'.
}
\description{
Load a matrix of OTUs in Qiime's format
}
\examples{
# see vignette
}
\seealso{
\code{\link{load_meta}} \code{\link{load_phenoData}}
}

