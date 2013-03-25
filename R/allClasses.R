setClass("MRexperiment", contains=c("eSet"), representation=representation(expSummary = "environment"),prototype = prototype( new( "VersionedBiobase",versions = c(classVersion("eSet"),MRexperiment = "1.0.0" ))))
            
setMethod("[", "MRexperiment", function (x, i, j, ..., drop = FALSE) {
        obj= callNextMethod()
        if(!missing(j)){
            obj@expSummary = new("environment",expSummary=as(expSummary(x)[j,1:2,...,drop=drop],"AnnotatedDataFrame"),cumNormStat=x@expSummary$cumNormStat)
            for(i in 1:length(pData(obj))){
                pData(obj)[,i] = factor(pData(obj)[,i])
            }
        }
        obj
})

newMRexperiment <- function(counts, phenoData=NULL, featureData=NULL,libSize=NULL, normFactors=NULL) {
    counts= as.matrix(counts)

    if( is.null( featureData ) )
      featureData <- annotatedDataFrameFrom(counts, byrow=TRUE)
    if( is.null( phenoData ) )
      phenoData   <- annotatedDataFrameFrom(counts, byrow=FALSE)
    if( is.null( libSize ) )
      libSize <- as.matrix(colSums(counts))
      rownames(libSize) = colnames(counts)
    if( is.null( normFactors ) ){
      normFactors <- as.matrix(rep( NA_real_, length(libSize) ))
      rownames(normFactors) = rownames(libSize)
    }

    obj <-new("MRexperiment", assayData = assayDataNew("environment",counts=counts),phenoData = phenoData,featureData = featureData ,expSummary = new("environment",expSummary=annotatedDataFrameFrom(counts,byrow=FALSE),cumNormStat=NULL))
    obj@expSummary$expSummary$libSize = libSize;
    obj@expSummary$expSummary$normFactors=normFactors;
        
    validObject( obj )
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

MRcounts <- function( obj ,norm=FALSE) {
   stopifnot( is( obj, "MRexperiment" ) )
   if(!norm){
    return(assayData(obj)[["counts"]])
   }
   if(any(is.na(normFactors(obj)))){
    return("Calculate the normalization factors first!")
   } else{
    sweep(assayData(obj)[["counts"]],2,as.vector(unlist(normFactors(obj)))/1000,"/")
   }
}

posterior.probs <- function( obj ) {
   stopifnot( is( obj, "MRexperiment" ) )
   assayData(obj)[["z"]]
}   

normFactors <- function( obj ) {
   stopifnot( is( obj, "MRexperiment" ) )
   nf <- pData(obj@expSummary$expSummary)[["normFactors"]]
   nf
}   

libSize<-function(obj){
   stopifnot( is( obj, "MRexperiment" ) )
   ls <- pData(obj@expSummary$expSummary)[["libSize"]]
   ls
}

expSummary<-function(obj){
  stopifnot( is( obj, "MRexperiment" ) )
  pData(obj@expSummary$expSummary)
}