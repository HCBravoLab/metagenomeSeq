#' Estimated effective samples per feature
#' 
#' Calculates the number of estimated effective samples per feature from the output
#' of a fitZig run. The estimated effective samples per feature is calculated as the
#' sum_1^n (n = number of samples) 1-z_i where z_i is the posterior probability a feature
#' belongs to the technical distribution.
#' 
#' @param obj The output of fitZig run on a MRexperiment object.
#' @return A list of the estimated effective samples per feature.
#' @seealso \code{\link{fitZig}} \code{\link{MRcoefs}} \code{\link{MRfulltable}}
#' 
calculateEffectiveSamples<-function(obj){
    rowSums(1-obj$z)
}
