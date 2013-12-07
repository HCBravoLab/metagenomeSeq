\name{filterData}
\alias{filterData}
\title{Filter datasets according to no. features present in features with at least a certain depth.}
\usage{
filterData(obj, present = 1, depth = 1000)
}
\arguments{
  \item{obj}{A MRexperiment object.}

  \item{present}{Features with at least 'present' postive
  samples.}

  \item{depth}{Sampls with at least this much depth of
  coverage}
}
\value{
A MRexperiment object.
}
\description{
Filter the data based on the number of present features
after filtering samples by depth of coverage. There are
many ways to filter the object, this is just one way.
}
\examples{
data(mouseData)
filterData(mouseData)
}

