\name{correlationTest}
\alias{corTest}
\alias{correlationTest}
\title{Pairwise correlation for each row of a matrix or MRexperiment object}
\usage{
correlationTest(obj, method = "pearson", alternative = "two.sided",
  norm = TRUE, log = TRUE, parallel = FALSE, cores = 2,
  override = FALSE, ...)
}
\arguments{
  \item{obj}{A MRexperiment object or count matrix.}

  \item{method}{One of 'pearson','spearman', or 'kendall'.}

  \item{alternative}{Indicates the alternative hypothesis
  and must be one of ‘two.sided’, ‘greater’ (positive) or
  ‘less’(negative). You can specify just the initial
  letter.}

  \item{norm}{Whether to aggregate normalized counts or not
  - if MRexperiment object.}

  \item{log}{Whether or not to log2 transform the counts -
  if MRexperiment object.}

  \item{parallel}{Parallelize the correlation testing?}

  \item{cores}{Number of cores to make use of if parallel
  == TRUE.}

  \item{override}{If the number of rows to test is over a
  thousand the test will not commence (unless
  override==TRUE).}

  \item{...}{Extra settings for mclapply.}
}
\value{
A matrix of size choose(number of rows, 2) by 2. The first
column corresponds to the correlation value. The second
column the p-value.
}
\description{
Calculates the pairwise correlation statistics and
associated p-values.
}
\examples{
data(mouseData)
cors = correlationTest(mouseData[1:10,],norm=FALSE,log=FALSE)
head(cors)
}

