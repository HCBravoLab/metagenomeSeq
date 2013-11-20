\name{plotGenus}
\alias{genusPlot}
\alias{plotGenus}
\title{Basic plot function of the raw or normalized data.}
\usage{
  plotGenus(obj, otuIndex, classIndex, log = TRUE,
    norm = TRUE, no = 1:length(otuIndex), labs = TRUE,
    xlab = NULL, ylab = NULL, jitter = TRUE,
    jitter.factor = 1, pch = 21, ret = FALSE, ...)
}
\arguments{
  \item{obj}{An MRexperiment object with count data.}

  \item{otuIndex}{A list of the otus with the same
  annotation.}

  \item{classIndex}{A list of the samples in their
  respective groups.}

  \item{log}{Whether or not to log2 transform the counts -
  if MRexperiment object.}

  \item{norm}{Whether or not to normalize the counts - if
  MRexperiment object.}

  \item{no}{Which of the otuIndex to plot.}

  \item{jitter.factor}{Factor value for jitter}

  \item{pch}{Standard pch value for the plot command.}

  \item{labs}{Whether to include group labels or not.
  (TRUE/FALSE)}

  \item{xlab}{xlabel for the plot.}

  \item{ylab}{ylabel for the plot.}

  \item{jitter}{Boolean to jitter the count data or not.}

  \item{ret}{Boolean to return the observed data that would
  have been plotted.}

  \item{...}{Additional plot arguments.}
}
\value{
  NA
}
\description{
  This function plots the abundance of a particular OTU by
  class. The function uses the estimated posterior
  probabilities to make technical zeros transparent.
}
\examples{
data(mouseData)
classIndex=list(controls=which(pData(mouseData)$diet=="BK"))
classIndex$cases=which(pData(mouseData)$diet=="Western")
otuIndex = grep("Strep",fData(mouseData)$fdata)
otuIndex=otuIndex[order(rowSums(MRcounts(mouseData)[otuIndex,]),decreasing=TRUE)]
plotGenus(mouseData,otuIndex,classIndex,no=1:2,xaxt="n",norm=FALSE,ylab="Strep normalized log(cpt)")
}
\seealso{
  \code{\link{cumNorm}}
}

