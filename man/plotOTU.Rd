\name{plotOTU}
\alias{plotOTU}
\title{Basic plot function of the raw or normalized data.}
\usage{
  plotOTU(obj, otu, classIndex, norm = TRUE,
    factor = 1, pch = 21, jitter = TRUE, ret = FALSE, ...)
}
\arguments{
  \item{obj}{A MRexperiment object with count data.}

  \item{otu}{The row number/OTU to plot.}

  \item{classIndex}{A list of the samples in their
  respective groups.}

  \item{norm}{Whether or not to normalize the counts.}

  \item{factor}{Factor value for jitter.}

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
# you can specify whether or not to normalize, and to what level
plotOTU(mouseData,otu=9083,classIndex,xlab="OTU log-normalized counts",norm=FALSE,xaxt="n",main="9083 feature abundances")
lablist<- c("Controls","Cases")
axis(1, at=seq(1,2,by=1), labels = lablist)
}
\seealso{
  \code{\link{cumNorm}}
}

