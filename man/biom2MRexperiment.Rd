\name{biom2MRexperiment}
\alias{biom2MRexperiment}
\title{Biome to MRexperiment objects}
\usage{
  biom2MRexperiment(obj)
}
\arguments{
  \item{obj}{The biome object file.}
}
\value{
  A MRexperiment object.
}
\description{
  Wrapper to convert biome files to MRexperiment objects.
}
\examples{
#library(biom)
#rich_dense_file = system.file("extdata", "rich_dense_otu_table.biom", package = "biom")
#x = read_biom(rich_dense_file)
#biom2MRexperiment(x)
}
\seealso{
  \code{\link{load_meta}} \code{\link{load_phenoData}}
  \code{\link{newMRexperiment}} \code{\link{load_biom}}
}

