#' Computes the weighted fold-change estimates and t-statistics.
#'
#' Wrapper to actually run the Expectation-maximization algorithm and estimate $f_{count}$ fits.
#' Maximum-likelihood estimates are approximated using the EM algorithm where we treat mixture membership $\deta_{ij}$ = 1 if $y_{ij}$
#' is generated from the zero point mass as latent indicator variables. The density is defined as $f_zig(y_{ij} = \pi_j(S_j) \cdot f_{0}(y_{ij}) 
#' +(1-\pi_j (S_j))\cdot f_{count}(y_{ij};\mu_i,\sigma_i^2)$.
#' The log-likelihood in this extended model is
#' $(1−\delta_{ij}) \log f_{count}(y;\mu_i,\sigma_i^2 )+\delta_{ij} \log \pi_j(s_j)+(1−\delta_{ij})\log (1−\pi_j (sj))$.
#' The responsibilities are defined as $z_{ij} = pr(\delta_{ij}=1 | data)$.
#'
#' @param obj An eSet object with count data.
#' @param mod The model for the count distribution.
#' @param s95 A vector of size M of the scaling values to be included in the model.
#' @param zeroMod The zero model, the model to account for the change in the number of OTUs observed as a linear effect of the depth of coverage.
#' @param useS95offset Boolean, whether to include the default scaling parameters in the model or not.
#' @param control The settings for fitZig. 
#' @param s The raw total counts for the various samples.
#' @return The fits, posterior probabilities, posterior probabilities used at time of convergence for each feature, ebayes (limma object) fit, among other data. 
#'
#' @name fitZig
#' @seealso \code{\link{cumNorm}} \code{\link{zigControl}}
#' @examples
#' model = model.matrix(~1+type+log2(s95/1000+1))
#' res = fitZig(obj = obj,mod=mod,useS95offset=FALSE)
fitZig <-
function(obj,mod,zeroMod=NULL,useS95offset=TRUE,control=zigControl()){

# Initialization
	tol = control$tol;
	maxit     = control$maxit;
	verbose   = control$verbose;
	
    stopifnot( is( obj, "MRexperiment" ) )
    if(any(is.na(normFactors(obj)))){return("At least one NA normalization factors")}
	
	y = MRcounts(obj)
	nc = ncol(y) #nsamples
	nr = nrow(y) #nfeatures

	zeroIndices=(y==0)
	z=matrix(0,nrow=nr, ncol=nc)
	z[zeroIndices]=0.5
	zUsed = z;
	curIt=0
	nllOld=rep(Inf, nr)
	nll=rep(Inf, nr)
	nllUSED=nll
	stillActive=rep(TRUE, nr);
	stillActiveNLL=rep(1, nr)
	
# Normalization step
		Nmatrix = log2(y+1)
		
		
# Initializing the model matrix
	if(useS95offset==TRUE){
		mmCount=cbind(mod,log2(as.matrix(normFactors(obj))/1000 +1))}
	else{ 
        mmCount=mod
    }

	if(is.null(zeroMod)){
        mmZero=model.matrix(~1+log(libSize(obj)))
    } else{ 
        mmZero=zeroMod 
    }
	
	modRank=ncol(mmCount);
# E-M Algorithm
		while(any(stillActive) && curIt<maxit) {
	
# M-step for count density (each feature independently)
			if(curIt==0){
				fit=doCountMStep(z, Nmatrix, mmCount, stillActive);
			} else {
				fit=doCountMStep(z, Nmatrix, mmCount, stillActive,fit)
			}

# M-step for zero density (all features together)
			zeroCoef = doZeroMStep(z, zeroIndices, mmZero)
			
# E-step
			z = doEStep(fit$residuals, zeroCoef$residuals, zeroIndices)
			zzdata<-getZ(z,zUsed,stillActive,nll,nllUSED);
			zUsed = zzdata$zUsed;
# NLL 
			nll = getNegativeLogLikelihoods(z, fit$residuals, zeroCoef$residuals)
			eps = getEpsilon(nll, nllOld)
			active = isItStillActive(eps, tol,stillActive,stillActiveNLL,nll)
			stillActive = active$stillActive;
			stillActiveNLL = active$stillActiveNLL;
			if(verbose==TRUE){
				cat(sprintf("it=%2d, nll=%0.2f, log10(eps+1)=%0.2f, stillActive=%d\n", curIt, mean(nll,na.rm=TRUE), log10(max(eps,na.rm=TRUE)+1), sum(stillActive)))
			}
			nllOld=nll
			curIt=curIt+1
    
			if(sum(rowSums((1-z)>0)<=modRank,na.rm=TRUE)>0){
				k = which(rowSums((1-z)>0)<=modRank)
				stillActive[k] = FALSE;
				stillActiveNLL[k] = nll[k]
			}
		}
	
        assayData(obj)[["z"]] <- z
        assayData(obj)[["zUsed"]] <- zUsed

		eb=limma::ebayes(fit$fit)
		dat = list(fit=fit$fit,countResiduals=fit$residuals,
				   z=z,eb=eb,taxa=rownames(obj),counts=y,zeroMod =mmZero,stillActive=stillActive,stillActiveNLL=stillActiveNLL,zeroCoef=zeroCoef)
		return(dat)
	}