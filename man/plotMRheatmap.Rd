\name{plotMRheatmap}
\alias{plotMRheatmap}
\title{Basic heatmap plot function for normalized counts.}
\usage{
  plotMRheatmap(obj, n, trials, log = TRUE, norm = TRUE,
    ...)
}
\arguments{
  \item{obj}{A MRexperiment object with count data.}

  \item{n}{The number of features to plot}

  \item{trials}{A vector of clinical information for.}

  \item{log}{Whether or not to log transform the counts.}

  \item{norm}{Whether or not to normalize the counts.}

  \item{...}{Additional plot arguments.}
}
\value{
  NA
}
\description{
  This function plots a heatmap of the "n" features with
  greatest variance across rows.
}
\examples{
data(mouseData)
trials = pData(mouseData)$diet
plotMRheatmap(obj=mouseData,n=200,trials=trials,cexRow = 0.4,cexCol = 0.4,trace="none")
}
\seealso{
  \code{\link{cumNormMat}}
}

