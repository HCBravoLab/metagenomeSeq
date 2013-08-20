\name{plotCorr}
\alias{plotCorr}
\title{Basic correlation plot function for normalized or unnormalized counts.}
\usage{
  plotCorr(obj, n, log = TRUE, norm = TRUE, fun = cor, ...)
}
\arguments{
  \item{obj}{A MRexperiment object with count data.}

  \item{n}{The number of features to plot}

  \item{log}{Whether or not to log transform the counts - if MRexperiment object.}

  \item{norm}{Whether or not to normalize the counts - if MRexperiment object.}

  \item{fun}{Function to calculate pair-wise relationships.
  Default is pearson correlation}

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
plotCorr(obj=mouseData,n=200,cexRow = 0.4,cexCol = 0.4,trace="none",dendrogram="none",
         col = colorRampPalette(brewer.pal(9, "RdBu"))(50))
}
\seealso{
  \code{\link{cumNormMat}}
}

