\name{aggregateByTaxonomy}
\alias{aggTax}
\alias{aggregateByTaxonomy}
\title{Aggregates a MRexperiment object by a particular taxonomic level.

Using the featureData information in the MRexperiment, calling aggregateByTaxonomy on a
MRexperiment and a particular featureData column (i.e. 'genus') will aggregate counts
to the desired level using the aggfun function (default colSums). Possible aggfun alternatives
include colMeans and colMedians.}
\usage{
aggregateByTaxonomy(obj, lvl, alternate = FALSE, norm = TRUE,
  aggfun = colSums)

aggTax(obj, lvl, alternate = FALSE, norm = TRUE, aggfun = colSums)
}
\arguments{
  \item{obj}{A MRexperiment object.}

  \item{lvl}{featureData column name from the MRexperiment
  object.}

  \item{alternate}{Use the rowname for undefined OTUs
  instead of aggregating to "no_match".}

  \item{norm}{Whether to aggregate normalized counts or
  not.}

  \item{aggfun}{Aggregation function.}
}
\value{
An aggregated count matrix.
}
\description{
Aggregates a MRexperiment object by a particular taxonomic
level.

Using the featureData information in the MRexperiment,
calling aggregateByTaxonomy on a MRexperiment and a
particular featureData column (i.e. 'genus') will aggregate
counts to the desired level using the aggfun function
(default colSums). Possible aggfun alternatives include
colMeans and colMedians.
}
\examples{
# not run
# aggregateByTaxonomy(mouseData,lvl="genus",norm=TRUE,aggfun=colMedians)
}

