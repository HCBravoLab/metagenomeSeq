\name{MRexperiment2biom}
\alias{MRexperiment2biom}
\title{MRexperiment to biom objects}
\usage{
MRexperiment2biom(obj, id = NULL)
}
\arguments{
  \item{obj}{The MRexperiment object.}

  \item{id}{Optional id for the biom matrix.}
}
\value{
A biom object.
}
\description{
Wrapper to convert MRexperiment objects to biom objects.
}
\seealso{
\code{\link{load_meta}} \code{\link{load_phenoData}}
\code{\link{newMRexperiment}} \code{\link{load_biom}}
\code{\link{biom2MRexperiment}}
}

