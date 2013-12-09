\name{MRcoefs}
\alias{MRcoefs}
\title{Table of top-ranked microbial marker gene from linear model fit}
\usage{
MRcoefs(obj, by = 2, coef = NULL, number = 10, taxa = obj$taxa,
  uniqueNames = FALSE, adjust.method = "fdr", group = 0, eff = 0,
  output = NULL)
}
\arguments{
  \item{obj}{A list containing the linear model fit
  produced by lmFit through fitZig.}

  \item{by}{Column number or column name specifying which
  coefficient or contrast of the linear model is of
  interest.}

  \item{coef}{Column number(s) or column name(s) specifying
  which coefficient or contrast of the linear model to
  display.}

  \item{number}{The number of bacterial features to pick
  out.}

  \item{taxa}{Taxa list.}

  \item{uniqueNames}{Number the various taxa.}

  \item{adjust.method}{Method to adjust p-values by.
  Default is "FDR". Options include "holm", "hochberg",
  "hommel", "bonferroni", "BH", "BY", "fdr", "none". See
  \code{\link{p.adjust}} for more details.}

  \item{group}{One of three choices, 0,1,2,3. 0: the sort
  is ordered by a decreasing absolute value coefficient
  fit. 1: the sort is ordered by the raw coefficient fit in
  decreasing order. 2: the sort is ordered by the raw
  coefficient fit in increasing order. 3: the sort is
  ordered by the p-value of the coefficient fit in
  increasing order.}

  \item{eff}{Restrict samples to have at least eff quantile
  effective samples.}

  \item{output}{Name of output file, including location, to
  save the table.}
}
\value{
Table of the top-ranked features determined by the linear
fit's coefficient.
}
\description{
Extract a table of the top-ranked features from a linear
model fit. This function will be updated soon to provide
better flexibility similar to limma's topTable.
}
\examples{
data(lungData)
k = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-k]
k = which(rowSums(MRcounts(lungTrim)>0)<10)
lungTrim = lungTrim[-k,]
cumNorm(lungTrim)
smokingStatus = pData(lungTrim)$SmokingStatus
mod = model.matrix(~smokingStatus)
settings = zigControl(maxit=1,verbose=FALSE)
fit = fitZig(obj = lungTrim,mod=mod,control=settings)
head(MRcoefs(fit))
}
\seealso{
\code{\link{fitZig}} \code{\link{MRtable}}
}

