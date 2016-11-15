#' Calculate the mixture proportions from the zero model / spike mass model
#' residuals.
#' 
#' F(x) = 1 / (1 + exp(-(x-m)/s)) (the CDF of the logistic distribution).
#' Provides the probability that a real-valued random variable X with a given
#' probability distribution will be found at a value less than or equal to x.
#' The output are the mixture proportions for the samples given the residuals
#' from the zero model.
#' 
#' 
#' @param residuals Residuals from the zero model.
#' @return Mixture proportions for each sample.
#' @seealso \code{\link{fitZig}}
getPi <- function(residuals) {
	plogis(residuals)
}
