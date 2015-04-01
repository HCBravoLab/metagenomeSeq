#' Depcrecated functions in the metagenomeSeq package.
#' 
#' These functions may be removed completely in the next release.
#' 
#' @usage deprecated_metagenomeSeq_function(x, value, ...)
#' @rdname metagenomeSeq-deprecated
#' @name metagenomeSeq-deprecated
#' @param x For assignment operators, the object that will undergo a replacement
#'  (object inside parenthesis).
#' @param value For assignment operators, the value to replace with 
#'  (the right side of the assignment).
#' @param ... For functions other than assignment operators, 
#'  parameters to be passed to the modern version of the function (see table).
#' @docType package
#' @export fitMeta
#' @aliases deprecated_metagenomeSeq_function fitMeta
#'
deprecated_metagenomeSeq_function <- function(x, value, ...){return(NULL)}
fitMeta <- function(...){.Deprecated("fitMeta",package="metagenomeSeq");return(fitLogNormal(...))}
