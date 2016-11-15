#' Computes the weighted fold-change estimates and t-statistics.
#' 
#' Wrapper to actually run the Expectation-maximization algorithm and estimate
#' $f_count$ fits.  Maximum-likelihood estimates are approximated using the EM
#' algorithm where we treat mixture membership $delta_ij = 1$ if $y_ij$ is
#' generated from the zero point mass as latent indicator variables. The
#' density is defined as $f_zig(y_ij = pi_j(S_j)*f_0(y_ij) +(1-pi_j (S_j)) *
#' f_count(y_ij; mu_i, sigma_i^2)$. The log-likelihood in this extended model
#' is: $(1-delta_ij) log f_count(y;mu_i,sigma_i^2 )+delta_ij log
#' pi_j(s_j)+(1-delta_ij) log (1-pi_j (s_j))$. The responsibilities are defined
#' as $z_ij = pr(delta_ij=1 | data)$.
#' 
#' 
#' @param obj A MRexperiment object with count data.
#' @param mod The model for the count distribution.
#' @param zeroMod The zero model, the model to account for the change in the
#' number of OTUs observed as a linear effect of the depth of coverage.
#' @param useCSSoffset Boolean, whether to include the default scaling
#' parameters in the model or not.
#' @param control The settings for fitZig.
#' @param useMixedModel Estimate the correlation between duplicate 
#' features or replicates using duplicateCorrelation.
#' @param ... Additional parameters for duplicateCorrelation.
#' @return A list of objects including:
#' \itemize{
#' 	\item{call - the call made to fitZig}
#' 	\item{fit  - 'MLArrayLM' Limma object of the weighted fit}
#' 	\item{countResiduals - standardized residuals of the fit}
#' 	\item{z - matrix of the posterior probabilities}
#' 	\item{eb - output of eBayes, moderated t-statistics, moderated F-statistics, etc}
#' 	\item{taxa - vector of the taxa names}
#' 	\item{counts - the original count matrix input}
#' 	\item{zeroMod - the zero model matrix}
#' 	\item{zeroCoef - the zero model fitted results}
#' 	\item{stillActive - convergence}
#' 	\item{stillActiveNLL - nll at convergence}
#' 	\item{dupcor - correlation of duplicates}
#' }
#' @export
#' @seealso \code{\link{cumNorm}} \code{\link{zigControl}}
#' @examples
#'
#' # This is a simple demonstration 
#' data(lungData)
#' k = grep("Extraction.Control",pData(lungData)$SampleType)
#' lungTrim = lungData[,-k]
#' k = which(rowSums(MRcounts(lungTrim)>0)<30)
#' lungTrim = cumNorm(lungTrim)
#' lungTrim = lungTrim[-k,]
#' smokingStatus = pData(lungTrim)$SmokingStatus
#' mod = model.matrix(~smokingStatus)
#' # The maxit is not meant to be 1 - this is for demonstration/speed
#' settings = zigControl(maxit=1,verbose=FALSE)
#' fit = fitZig(obj = lungTrim,mod=mod,control=settings)
#' 
fitZig <- function(obj, 
                   mod, 
                   zeroMod=NULL,
                   useCSSoffset=TRUE,
                   control=zigControl(),
                   useMixedModel=FALSE,
                   ...)
{
	stopifnot( is( obj, "MRexperiment" ) )
	if(any(is.na(normFactors(obj)))) stop("At least one NA normalization factors")
	if(any(is.na(libSize(obj)))) stop("Calculate the library size first!")
	
	y <- MRcounts(obj, norm=FALSE, log=FALSE)
	nc <- ncol(y) #nsamples
	nr <- nrow(y) #nfeatures

  # Normalization step
	Nmatrix <- log2(y + 1)
		
  # Initializing the model matrix
	if (useCSSoffset == TRUE){
		if (any(is.na(normFactors(obj)))) {
		  stop("Calculate the normalization factors first!")
		}
		mmCount <- cbind(mod, log2(normFactors(obj)/1000 + 1))
		colnames(mmCount)[ncol(mmCount)] <- "scalingFactor"
	} else { 
		mmCount <- mod
	}

	if (is.null(zeroMod)) {
	  if (any(is.na(libSize(obj)))) {
		  stop("Calculate the library size first!") 
		}
		
	  mmZero <- model.matrix(~1+log(libSize(obj)))
	} else { 
		mmZero <- zeroMod
	}
	
	dat <- .do_fitZig(Nmatrix, mmCount, mmZero, control=control, useMixedModel=useMixedModel, ...)

	assayData(obj)[["z"]] <- dat$z
	assayData(obj)[["zUsed"]] <- dat$zUsed
	dat$zUsed <- NULL
	
	dat <- c(dat, list(call=match.call(),taxa=rownames(obj),counts=y))
	dat
}

.do_fitZig <- function(y, 
                       count_model_matrix, 
                       zero_model_matrix, 
                       control=zigControl(), 
                       useMixedModel=FALSE,
                       ...) 
{
  # Initialization
  tol <- control$tol
  maxit <- control$maxit
  verbose <- control$verbose
  dfMethod <- control$dfMethod
  pvalMethod <- control$pvalMethod
  
  nr <- nrow(y)
  nc <- ncol(y)
  
  zeroIndices <- (y == 0)
  z <- matrix(0, nrow=nr, ncol=nc)
  z[zeroIndices] <- 0.5
  zUsed <- z
  
  curIt <- 0
  nllOld <- rep(Inf, nr)
  nll <- rep(Inf, nr)
  nllUSED <- nll
  stillActive <- rep(TRUE, nr)
  stillActiveNLL <- rep(1, nr)
  dupcor <- NULL
  
  modRank <- ncol(count_model_matrix)
  # E-M Algorithm
  while (any(stillActive) && (curIt < maxit)) {
    
    # M-step for count density (each feature independently)
    if (curIt == 0) {
      fit <- doCountMStep(z, y, count_model_matrix, stillActive, dfMethod=dfMethod)
    } else {
      fit <- doCountMStep(z, y, count_model_matrix, stillActive, fit2=fit, dfMethod=dfMethod)
    }
    
    # M-step for zero density (all features together)
    zeroCoef <- doZeroMStep(z, zeroIndices, zero_model_matrix)
    
    # E-step
    z <- doEStep(fit$residuals, zeroCoef$residuals, zeroIndices)
    zzdata <- getZ(z, zUsed, stillActive, nll, nllUSED);
    zUsed <- zzdata$zUsed;
    
    # NLL 
    nll <- getNegativeLogLikelihoods(z, fit$residuals, zeroCoef$residuals)
    eps <- getEpsilon(nll, nllOld)
    active <- isItStillActive(eps, tol,stillActive,stillActiveNLL,nll)
    stillActive <- active$stillActive;
    stillActiveNLL <- active$stillActiveNLL;
    if (verbose == TRUE){
      cat(sprintf("it=%2d, nll=%0.2f, log10(eps+1)=%0.2f, stillActive=%d\n", curIt, mean(nll,na.rm=TRUE), log10(max(eps,na.rm=TRUE)+1), sum(stillActive)))
    }
    nllOld <- nll
    curIt <- curIt + 1
    
    if (sum(rowSums((1-z) > 0) <= modRank, na.rm=TRUE) > 0) {
      k <- which(rowSums((1-z) > 0) <= modRank)
      stillActive[k] <- FALSE;
      stillActiveNLL[k] <- nll[k]
    }
  }
  
  if (useMixedModel == TRUE) {
    dupcor <- duplicateCorrelation(y, count_model_matrix, weights=(1-z), ...)
    fit$fit <- limma::lmFit(y, count_model_matrix, weights=(1-z), correlation=dupcor$consensus, ...)
    countCoef <- fit$fit$coefficients
    countMu <- tcrossprod(countCoef, count_model_matrix)
    fit$residuals <- sweep((y-countMu), 1, fit$fit$sigma, "/")
  }
  
  eb <- limma::eBayes(fit$fit)
  dat <- list(fit=fit$fit, countResiduals=fit$residuals,
              z=z, zUsed=zUsed, eb=eb, zeroMod=zero_model_matrix, stillActive=stillActive, 
              stillActiveNLL=stillActiveNLL, zeroCoef=zeroCoef, dupcor=dupcor)
  dat
}

# #' Function to perform fitZig bootstrap
# #' 
# #' Calculates bootstrap stats
# #' 
# #' @param y Log-transformed matrix
# #' @param y string for the y-axis
# #' @param norm is the data normalized?
# #' @param log is the data logged?
# #' @return vector of x,y labels
# #' 
# performBoostrap<-function(fit){


# 	zeroIndices=(y==0)
# 	z=matrix(0,nrow=nr, ncol=nc)
# 	z[zeroIndices]=0.5
# 	zUsed = z

# 	curIt=0
# 	nllOld=rep(Inf, nr)
# 	nll=rep(Inf, nr)
# 	nllUSED=nll
# 	stillActive=rep(TRUE, nr)
# 	stillActiveNLL=rep(1, nr)
	
# 	tt <- fit$fit$coef[,coef] / fit$fit$stdev.unscaled[,coef] / fit$fit$sigma
# 	perms = replicate(B,sample(mmCount[,coef]))
# 	mmCount1=mmCount[,-coef]

# # Normalization step
# 	Nmatrix = log2(y+1)
		
# # Initializing the model matrix
# 	if(useCSSoffset==TRUE){
# 		if(any(is.na(normFactors(obj)))){stop("Calculate the normalization factors first!")}
# 		mmCount=cbind(mod,log2(normFactors(obj)/1000 +1))
# 		colnames(mmCount)[ncol(mmCount)] = "scalingFactor"
# 	}
# 	else{ 
# 		mmCount=mod
# 	}

# 	if(is.null(zeroMod)){
# 		if(any(is.na(libSize(obj)))){ stop("Calculate the library size first!") }
# 			mmZero=model.matrix(~1+log(libSize(obj)))
# 		} else{ 
# 			mmZero=zeroMod
# 		}
	
# 	modRank=ncol(mmCount)
# 	# E-M Algorithm
# 	while(any(stillActive) && curIt<maxit) {
	
# # M-step for count density (each feature independently)
# 		if(curIt==0){
# 			fit=doCountMStep(z, Nmatrix, mmCount, stillActive,dfMethod=dfMethod);
# 		} else {
# 			fit=doCountMStep(z, Nmatrix, mmCount, stillActive,fit2=fit,dfMethod=dfMethod)
# 		}

# # M-step for zero density (all features together)
# 		zeroCoef = doZeroMStep(z, zeroIndices, mmZero)
			
# # E-step
# 		z = doEStep(fit$residuals, zeroCoef$residuals, zeroIndices)
# 		zzdata<-getZ(z,zUsed,stillActive,nll,nllUSED);
# 		zUsed = zzdata$zUsed;
# # NLL 
# 		nll = getNegativeLogLikelihoods(z, fit$residuals, zeroCoef$residuals)
# 		eps = getEpsilon(nll, nllOld)
# 		active = isItStillActive(eps, tol,stillActive,stillActiveNLL,nll)
# 		stillActive = active$stillActive;
# 		stillActiveNLL = active$stillActiveNLL;
# 		if(verbose==TRUE){
# 			cat(sprintf("it=%2d, nll=%0.2f, log10(eps+1)=%0.2f, stillActive=%d\n", curIt, mean(nll,na.rm=TRUE), log10(max(eps,na.rm=TRUE)+1), sum(stillActive)))
# 		}
# 		nllOld=nll
# 		curIt=curIt+1

# 		if(sum(rowSums((1-z)>0)<=modRank,na.rm=TRUE)>0){
# 			k = which(rowSums((1-z)>0)<=modRank)
# 			stillActive[k] = FALSE;
# 			stillActiveNLL[k] = nll[k]
# 		}
# 	}
# }
