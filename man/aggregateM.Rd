\name{aggregateM}
\alias{aggregateM}
\title{Aggregates counts by a particular classification.}
\usage{
  aggregateM(obj, taxa, lvl, split = ";")
}
\arguments{
  \item{obj}{A MRexperiment object.}

  \item{lvl}{The level to go up (numeric, 1,2,3).}

  \item{taxa}{A vector of taxa annotations with splits}

  \item{split}{The way character strings in taxa in the obj
  are split.}
}
\value{
  Updated object with counts aggregated to the various
  taxanomic levels.
}
\description{
  This function takes a MRexperiment object of data at a
  particular level with feature information allowing for
  aggregation of counts to a particular level. This method
  assumes taxa begin at the highest level and continue to
  the current level, reverse assumes taxa begin at the
  lowest level.
}

