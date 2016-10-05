setClass("MRexperiment", contains=c("eSet"), representation=representation(expSummary = "environment"),prototype = prototype( new( "VersionedBiobase",versions = c(classVersion("eSet"),MRexperiment = "1.0.0" ))))
            
setMethod("[", "MRexperiment", function (x, i, j, ..., drop = FALSE) {
        obj= callNextMethod()
        if(!missing(j)){
            obj@expSummary = new("environment",expSummary=as(expSummary(x)[j,1:2,...,drop=drop],"AnnotatedDataFrame"),cumNormStat=x@expSummary$cumNormStat)
            if(length(pData(obj))>0){
              for(i in 1:length(pData(obj))){
                if(is.factor(pData(obj)[,i])){
                  pData(obj)[,i] = factor(pData(obj)[,i])
                } else {
                  pData(obj)[,i] = pData(obj)[,i]
                }
              }
            }
        }
        obj
})
setMethod("colSums", signature ="MRexperiment", function (x, ...) {
    callNextMethod(MRcounts(x),...)
})
setMethod("rowSums", signature="MRexperiment", function (x, ...) {
    callNextMethod(MRcounts(x),...)
})
setMethod("rowMeans", signature="MRexperiment", function (x, ...) {
    callNextMethod(MRcounts(x),...)
})
setMethod("colMeans", signature="MRexperiment", function (x, ...) {
    callNextMethod(MRcounts(x),...)
})

#' Access the normalization factors in a MRexperiment object
#'
#' Function to access/replace the scaling factors, aka the normalization factors, of
#' samples in a MRexperiment object.
#'
#' @name normFactors
#' @aliases normFactors,MRexperiment-method normFactors normFactors<-
#' @docType methods
#' @param obj a \code{MRexperiment} object.
#' @return Normalization scaling factors
#' @author Joseph N. Paulson
#' @examples
#'
#' data(lungData)
#' head(normFactors(lungData))
#'
setGeneric("normFactors",function(object){standardGeneric("normFactors")}) 
setGeneric("normFactors<-",function(object,value){standardGeneric("normFactors<-")})

setMethod("normFactors", signature="MRexperiment",function(object) {
   nf <- expSummary(object)$normFactors
   nf <- as.data.frame(nf)
   colnames(nf) <- "normFactors"
   rownames(nf) <- sampleNames(object)
   nf
 })

setReplaceMethod("normFactors", signature=c(object="MRexperiment", value="numeric"),
  function( object, value ) {
   pData(object@expSummary$expSummary)$normFactors <- value
   validObject( object )
   object
})

#' Access sample depth of coverage from MRexperiment object
#'
#' Access or replace the libSize vector represents the column (sample specific) sums of features,
#' i.e. the total number of reads for a sample or depth of coverage. It is used by
#' \code{\link{fitZig}}.
#'
#' @name libSize
#' @aliases libSize,MRexperiment-method libSize libSize<-
#' @docType methods
#' @param obj a \code{MRexperiment} object.
#' @return Library sizes
#' @author Joseph N. Paulson, jpaulson@@umiacs.umd.edu
#' @examples
#'
#' data(lungData)
#' head(libSize(lungData))
#'
setGeneric("libSize",function(object){standardGeneric("libSize")})
setGeneric("libSize<-",function(object,value){standardGeneric("libSize<-")})

setMethod("libSize", signature="MRexperiment",function(object) {
   ls <- expSummary(object)$libSize
   ls <- as.data.frame(ls)
   colnames(ls) <- "libSize"
   rownames(ls) <- sampleNames(object)
   ls
 })

setReplaceMethod("libSize", signature=c(object="MRexperiment", value="numeric"),
  function( object, value ) {
   pData(object@expSummary$expSummary)$libSize <- value
   validObject( object )
   object
})

#' Create a MRexperiment object
#' 
#' This function creates a MRexperiment object from a matrix or data frame of
#' count data.
#' 
#' See \code{\link{MRexperiment-class}} and \code{eSet} (from the Biobase
#' package) for the meaning of the various slots.
#' 
#' @param counts A matrix or data frame of count data. The count data is
#' representative of the number of reads annotated for a feature (be it gene,
#' OTU, species, etc). Rows should correspond to features and columns to
#' samples.
#' @param phenoData An AnnotatedDataFrame with pertinent sample information.
#' @param featureData An AnnotatedDataFrame with pertinent feature information.
#' @param libSize libSize, library size, is the total number of reads for a
#' particular sample.
#' @param normFactors normFactors, the normalization factors used in either the
#' model or as scaling factors of sample counts for each particular sample.
#' @return an object of class MRexperiment
#' @author Joseph N Paulson
#' @examples
#' 
#' cnts = matrix(abs(rnorm(1000)),nc=10)
#' obj <- newMRexperiment(cnts)
#' 
newMRexperiment <- function(counts, phenoData=NULL, featureData=NULL,libSize=NULL, normFactors=NULL) {
    counts= as.matrix(counts)

    if( is.null( featureData ) ){
      featureData <- annotatedDataFrameFrom(counts, byrow=TRUE)
    }
    if( is.null( phenoData ) ){
      phenoData   <- annotatedDataFrameFrom(counts, byrow=FALSE)
    }
    if( is.null( libSize ) ){
      libSize <- as.matrix(colSums(counts))
      rownames(libSize) = colnames(counts)
    }
    if( is.null( normFactors ) ){
      normFactors <- as.matrix(rep( NA_real_, length(libSize) ))
      rownames(normFactors) = rownames(libSize)
    }

    obj <-new("MRexperiment", assayData = assayDataNew("environment",counts=counts),phenoData = phenoData,featureData = featureData ,expSummary = new("environment",expSummary=annotatedDataFrameFrom(counts,byrow=FALSE),cumNormStat=NULL))
    obj@expSummary$expSummary$libSize = libSize;
    obj@expSummary$expSummary$normFactors=normFactors;
    validObject(obj)
    obj
}
setValidity( "MRexperiment", function( object ) {
    if( is.null(assayData(object)$counts))
        return( "There are no counts!" )
#    if( ncol(MRcounts(object)) != length(normFactors(object)))
#        return( "Experiment summary got hacked!" )
#    if( ncol(MRcounts(object)) != length(libSize(object)))
#        return( "Experiment summary got hacked!" )
    TRUE
} )
#' Accessor for the counts slot of a MRexperiment object
#' 
#' The counts slot holds the raw count data representing (along the rows) the
#' number of reads annotated for a particular feature and (along the columns)
#' the sample.
#' 
#' 
#' @name MRcounts
#' @aliases MRcounts,MRexperiment-method MRcounts
#' @docType methods
#' @param obj a \code{MRexperiment} object.
#' @param norm logical indicating whether or not to return normalized counts.
#' @param log TRUE/FALSE whether or not to log2 transform scale.
#' @param sl The value to scale by (default=1000).
#' @return Normalized or raw counts
#' @author Joseph N. Paulson, jpaulson@@umiacs.umd.edu
#' @examples
#' 
#' data(lungData)
#' head(MRcounts(lungData))
#' 
MRcounts <- function(obj,norm=FALSE,log=FALSE,sl=1000) {
   stopifnot( is( obj, "MRexperiment" ) )
   if(!norm){
    x=assayData(obj)[["counts"]]
   }
   else{
    if(any(is.na(normFactors(obj)))){
      x=cumNormMat(obj,sl=sl)
    } else{
      x=sweep(assayData(obj)[["counts"]],2,as.vector(unlist(normFactors(obj)))/sl,"/")
    }
   }
   if(!log){
    return(x)
   } else{
    return(log2(x+1))
   }
}
#' Access the posterior probabilities that results from analysis
#' 
#' Accessing the posterior probabilities following a run through
#' \code{\link{fitZig}}
#' 
#' 
#' @name posteriorProbs
#' @aliases posteriorProbs,MRexperiment-method posteriorProbs
#' @docType methods
#' @param obj a \code{MRexperiment} object.
#' @return Matrix of posterior probabilities
#' @author Joseph N. Paulson, jpaulson@@umiacs.umd.edu
#' @examples
#' 
#' # see vignette
#' 
posteriorProbs <- function( obj ) {
   stopifnot( is( obj, "MRexperiment" ) )
   assayData(obj)[["z"]]
}

#' Access MRexperiment object experiment data
#' 
#' The expSummary vectors represent the column (sample specific) sums of
#' features, i.e. the total number of reads for a sample, libSize and also the
#' normalization factors, normFactor.
#' 
#' 
#' @name expSummary
#' @aliases expSummary,MRexperiment-method expSummary
#' @docType methods
#' @param obj a \code{MRexperiment} object.
#' @return Experiment summary table
#' @author Joseph N. Paulson, jpaulson@@umiacs.umd.edu
#' @examples
#' 
#' data(mouseData)
#' expSummary(mouseData)
#' 
expSummary<-function(obj){
  stopifnot( is( obj, "MRexperiment" ) )
  pData(obj@expSummary$expSummary)
}
#' Check if MRexperiment or matrix and return matrix
#'
#' Function to check if object is a MRexperiment
#' class or matrix 
#'
#' @name returnAppropriateObj
#' @param obj a \code{MRexperiment} or \code{matrix} object
#' @param norm return a normalized \code{MRexperiment} matrix
#' @param log return a log transformed \code{MRexperiment} matrix
#' @param sl scaling value
#' @return Matrix
#' @examples
#'
#' data(lungData)
#' head(returnAppropriateObj(lungData,norm=FALSE,log=FALSE))
#'
returnAppropriateObj <- function(obj,norm,log,sl=1000) {
  if(class(obj)=="MRexperiment"){
    mat = MRcounts(obj,norm=norm,log=log,sl=sl)
  } else if(class(obj) == "matrix") {
    mat = obj
  } else {
    stop("Object needs to be either a MRexperiment object or matrix")
  }
  mat
}
