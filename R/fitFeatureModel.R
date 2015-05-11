#' Computes differential abundance analysis using a zero-inflated log-normal model
#' 
#' Wrapper to actually run zero-inflated log-normal model given a MRexperiment object
#' and model matrix. User can decide to shrink parameter estimates.
#' 
#' @param obj A MRexperiment object with count data.
#' @param mod The model for the count distribution.
#' @param coef Coefficient of interest to grab log fold-changes.
#' @param B Number of iterations to perform. If >1, performs permutation test.
#' @param szero TRUE/FALSE, shrink zero component parameters.
#' @param spos TRUE/FALSE, shrink positive component parameters.
#' @return A list of objects including:
#' \itemize{
#'  \item{call - the call made to fitFeatureModel}
#'  \item{fitZeroLogNormal  - list of parameter estimates for the zero-inflated log normal model}
#'  \item{pvalues - calculated p-values}
#'  \item{permuttedfits - permutted z-score estimates under the null}
#' }
#' @seealso \code{\link{cumNorm}}
#' @examples
#' 
#' data(lungData)
#' lungData = lungData[,-which(is.na(pData(lungData)$SmokingStatus))]
#' lungData=filterData(lungData,present=5,depth=1)
#' lungData <- cumNorm(lungData, p=.5)
#' s <- normFactors(lungData)
#' pd <- pData(lungData)
#' pd <- cbind(pd,norm=log(s/median(s)))
#' mod <- model.matrix(~1+SmokingStatus, data=pd)
#' lungres1 = fitFeatureModel(lungData,mod)
#' 
fitFeatureModel<-function(obj,mod,coef=2,B=1,szero=FALSE,spos=TRUE){

  stopifnot(is(obj, "MRexperiment"))
  if (any(is.na(normFactors(obj)))) 
      stop("At least one NA normalization factors")
  if (any(is.na(libSize(obj)))) 
      stop("Calculate the library size first!")    
  if (any(is.na(normFactors(obj)))) {
      stop("Calculate the normalization factors first!")
  }
  mmCount = cbind(mod, log(normFactors(obj)/median(normFactors(obj))))
  colnames(mmCount)[ncol(mmCount)] = "scalingFactor"
  
  if(ncol(mmCount)>3){ stop("Can't analyze currently.") }

  # These pieces get to be a part of the new zero-ln model!
  show("Running zero-ln model")
  fitzeroln = fitZeroLogNormal(obj,mmCount,coef=coef,szero=szero,spos=spos)
  zscore = fitzeroln$logFC/fitzeroln$se

  if(B>1){
    show("Running permutations")
    permutations = replicate(B,sample(mmCount[,coef]))
    mmCountPerm  = mmCount
    export=c("fitZeroLogNormal","calcPosComponent",
      "calcZeroComponent","calcShrinkParameters",
      "calcZeroAdjustment","calcStandardError")
    
    permuttedFits = foreach(i = seq(B),.export=export,.errorhandling="remove",
      .packages=c("metagenomeSeq","glmnet")) %dopar% {
        mmCountPerm[,coef] = permutations[,i]
        permFit = fitZeroLogNormal(obj,mmCountPerm,coef=coef,szero=szero,spos=spos)
        permFit$logFC/permFit$se
      }

    zperm = abs(sapply(permuttedFits,function(i)i))
    pvals = rowMeans(zperm>=abs(zscore),na.rm=TRUE)
  } else {
    permuttedFits = NULL
    pvals = 2*(1-pnorm(abs(zscore)))
  }
  res = list(call=match.call(),fitZeroLogNormal=fitzeroln,pvalues=pvals,permuttedFits=permuttedFits)
  res
}