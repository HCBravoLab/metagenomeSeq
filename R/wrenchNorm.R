#' Computes normalization factors using wrench instead of cumNorm
#' 
#' Behaves in a similar manner to cumNorm but uses method published by
#' M. Sentil Kumar et al. (2018) to compute normalization factors which consider
#' compositional bias introduced by sequencers
#' 
#' @param obj an MRexperiment object
#' @param condition case control label that wrench uses to calculate normalization factors
#' @seealso \code{\link{cumNorm}}
#'
wrenchNorm <- function(obj, condition) {
  count_data <- MRcounts(obj, norm = FALSE)
  W <- wrench(count_data, condition = condition)
  normFactors(obj) <- W$nf
  return(obj)
}