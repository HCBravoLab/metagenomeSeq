\name{plotBubble}
\alias{plotBubble}
\title{Basic plot of binned vectors.}
\usage{
plotBubble(yvector, xvector, sigvector = NULL, nbreaks = 10, ret = FALSE,
  scale = 1, ...)
}
\arguments{
  \item{yvector}{A vector of values represented along
  y-axis.}

  \item{xvector}{A vector of values represented along
  x-axis.}

  \item{sigvector}{A vector of the names of significant
  features (names should match x/yvector).}

  \item{nbreaks}{Number of bins to break yvector and
  xvector into.}

  \item{ret}{Boolean to return the observed data that would
  have been plotted.}

  \item{scale}{Scaling of circle bins.}

  \item{...}{Additional plot arguments.}
}
\value{
NA
}
\description{
This function plots takes two vectors, calculates the
contingency table and plots circles sized by the
contingency table value. Optional significance vectors of
the values significant will shade the circles by proportion
of significance.
}
\examples{
data(mouseData)
mouseData = mouseData[which(rowSums(mouseData)>139),]
sparsity = rowMeans(MRcounts(mouseData)==0)
lor = log(fitPA(mouseData,cl=pData(mouseData)[,3])$oddsRatio)
plotBubble(lor,sparsity,main="lor ~ sparsity")
}
\seealso{
\code{\link{plotMRheatmap}}
}

