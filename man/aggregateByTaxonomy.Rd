\name{aggregateByTaxonomy}
\alias{aggTax}
\alias{aggregateByTaxonomy}
\alias{aggregateTax,}
\title{Aggregates a MRexperiment object by a particular taxonomic level.}
\usage{
aggregateByTaxonomy(obj, lvl, norm = TRUE, aggfun = colSums)
}
\arguments{
  \item{obj}{A MRexperiment object.}

  \item{lvl}{featureData column name from the MRexperiment
  object.}

  \item{norm}{Whether to aggregate normalized counts or
  not.}

  \item{aggfun}{Aggregation function.}
}
\value{
An aggregated count matrix.
}
\description{
Using the featureData information in the MRexperiment,
calling aggregateByTaxonomy on a MRexperiment and a
particular featureData column (i.e. 'genus') will aggregate
counts to the desired level by with the aggfun function
(default colSums). Possible aggfun alternatives include
colMeans and colMedians.
}
\examples{
# not run
# aggregateByTaxonomy(mouseData,lvl="genus",norm=TRUE,aggfun=colMedians)
}

