#' Computes normalization factors using wrench instead of cumNorm
#' 
#' Calculates normalization factors using method published by
#' M. Sentil Kumar et al. (2018) to compute normalization factors which 
#' considers compositional bias introduced by sequencers.
#' 
#' @param obj an MRexperiment object
#' @param condition case control label that wrench uses to calculate 
#' normalization factors
#' @return an MRexperiment object with updated normalization factors.
#' Accessible by \code{\link{normFactors}}.
#' @seealso \code{\link{cumNorm}} \code{\link{fitZig}}
#'
#' @examples
#' 
#' data(mouseData)
#' mouseData <- wrenchNorm(mouseData, condition = mouseData$diet)
#' head(normFactors(mouseData))
#' 
wrenchNorm <- function(obj, condition) {
  count_data <- MRcounts(obj, norm = FALSE)
  W <- wrench(count_data, condition = condition)
  normFactors(obj) <- W$nf
  return(obj)
}