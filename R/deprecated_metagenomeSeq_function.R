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
#' @aliases deprecated_metagenomeSeq_function fitMeta load_phenoData load_meta load_biom load_metaQ
#'
deprecated_metagenomeSeq_function <- function(x, value, ...){return(NULL)}
fitMeta <- function(...){.Deprecated("fitMeta",package="metagenomeSeq");return(fitLogNormal(...))}
load_phenoData <- function(...){.Deprecated("load_phenoData",package="metagenomeSeq");return(loadPhenodata(...))}
load_biom <- function(...){.Deprecated("load_biom",package="metagenomeSeq");return(loadBiom(...))}
load_meta <- function(...){.Deprecated("load_meta",package="metagenomeSeq");return(loadMeta(...))}
load_metaQ <- function(...){.Deprecated("load_metaQ",package="metagenomeSeq");return(loadMetaQ(...))}
