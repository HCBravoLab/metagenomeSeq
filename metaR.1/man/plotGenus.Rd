\name{plotGenus}
\alias{genusPlot}
\alias{plotGenus}
\title{Basic plot function of the raw or normalized data.}
\usage{
  plotGenus(obj, otuIndex, classIndex, norm = TRUE, no=1:length(otuIndex), jitter = TRUE, factor = 1,
    pch = 21, ret = FALSE, ...)
}
\arguments{
  \item{obj}{An MRexperiment object with count data.}

  \item{otuIndex}{A list of the otus with the same
  annotation.}

  \item{classIndex}{A list of the samples in their
  respective groups.}

  \item{norm}{Whether or not to normalize the counts.}

  \item{no}{Which of the otuIndex to plot.}

  \item{factor}{Factor value for jitter}

  \item{pch}{Standard pch value for the plot command.}

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
plotGenus(mouseData,otuIndex,classIndex,xlab="OTU log-normalized counts",no=1:2,xaxt="n",norm=FALSE,ylab="Strep")
lablist<-rep(c("Controls","Cases"),times=2)
axis(1, at=seq(1,4,by=1), labels = lablist)

}
\seealso{
  \code{\link{cumNorm}}
}

