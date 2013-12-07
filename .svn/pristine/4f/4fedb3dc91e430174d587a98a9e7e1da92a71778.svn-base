\name{plotOrd}
\alias{plotOrd}
\title{Plot of either PCA or MDS coordinates for the distances of normalized or unnormalized counts.}
\usage{
plotOrd(obj, tran = FALSE, comp = 1:2, log = TRUE, norm = TRUE,
  usePCA = TRUE, useDist = FALSE, distfun = stats::dist,
  dist.method = "euclidian", ret = FALSE, n = NULL, ...)
}
\arguments{
  \item{obj}{A MRexperiment object or count matrix.}

  \item{tran}{Transpose the matrix.}

  \item{comp}{Which components to display}

  \item{usePCA}{TRUE/FALSE whether to use PCA or MDS
  coordinates (TRUE is PCA).}

  \item{useDist}{TRUE/FALSE whether to calculate
  distances.}

  \item{distfun}{Distance function, default is stats::dist}

  \item{dist.method}{If useDist==TRUE, what method to
  calculate distances.}

  \item{log}{Whether or not to log2 the counts - if
  MRexperiment object.}

  \item{norm}{Whether or not to normalize the counts - if
  MRexperiment object.}

  \item{ret}{Whether or not to output the coordinates.}

  \item{n}{Number of features to make use of in calculating
  your distances.}

  \item{...}{Additional plot arguments.}
}
\value{
NA
}
\description{
This function plots the PCA / MDS coordinates for the "n"
features of interest. Potentially uncovering batch effects
or feature relationships.
}
\examples{
data(mouseData)
cl = pData(mouseData)[,3]
plotOrd(mouseData,tran=TRUE,useDist=TRUE,pch=21,bg=factor(cl),usePCA=FALSE)
}
\seealso{
\code{\link{cumNormMat}}
}

