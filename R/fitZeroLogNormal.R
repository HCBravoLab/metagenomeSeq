#' Compute the log fold-change estimates for the zero-inflated log-normal model
#' 
#' Run the zero-inflated log-normal model given a MRexperiment object
#' and model matrix. Not for the average user, assumes structure of the model matrix. 
#' 
#' @param obj A MRexperiment object with count data.
#' @param mod The model for the count distribution.
#' @param coef Coefficient of interest to grab log fold-changes.
#' @param szero TRUE/FALSE, shrink zero component parameters.
#' @param spos TRUE/FALSE, shrink positive component parameters.
#' @return A list of objects including:
#' \itemize{
#'  \item{logFC - the log fold-change estimates}
#'  \item{adjFactor  - the adjustment factor based on the zero component}
#'  \item{se - standard error estimates}
#'  \item{fitln - parameters from the log-normal fit}
#'  \item{fitzero - parameters from the logistic fit}
#'  \item{zeroRidge - output from the ridge regression}
#'  \item{posRidge - output from the ridge regression}
#'  \item{tauPos - estimated tau^2 for positive component}
#'  \item{tauZero - estimated tau^2 for zero component}
#'  \item{exclude - features to exclude for various reasons, e.g. all zeros}
#'  \item{zeroExclude - features to exclude for various reasons, e.g. all zeros}
#' }
#' @seealso \code{\link{cumNorm}} \code{\link{fitFeatureModel}}
fitZeroLogNormal<-function(obj,mod,coef=2,szero=TRUE,spos=TRUE){
  positiveMod = mod[,-ncol(mod)]
  zeroMod = mod
  
  nf <- normFactors(obj)
  mat <- MRcounts(obj, norm=TRUE, log=FALSE,sl=median(nf))
  posIndices = mat>0

  nr = nrow(mat)
  nc = ncol(mat)
  exclude = zeroExclude = tauZero = tauPos = posRidge = zeroRidge = NULL

  results = array(NA,dim=c(nr,3))
  rownames(results) = rownames(mat)
  colnames(results) = c("logFC","adjFactor","se")

  # calc log-normal component
  fitln = calcPosComponent(mat,positiveMod,posIndices)

  # Don't calculate shrinkage with special cases
  zeros2 = which(fitln[,"s2"]==0)
  rs = rowsum(t(1-(1-posIndices)),positiveMod[,coef])
  exclude = union(which(rs[1,]<=1),which(rs[2,]<=1))
  zeroExclude  = which(colSums(rs)>=(nc-3))
  exclude = union(zeros2,exclude); if(length(exclude)==0) exclude=NULL
  if(length(zeroExclude)==0) zeroExclude=NULL

  sdensity = density(fitln[,"s2"],na.rm=TRUE)
  smode = sdensity$x[which.max(sdensity$y)]
  if(length(zeros2)>0) fitln[zeros2,"s2"] = smode

  # shrink positive
  if(spos==TRUE){
    shrinkPos<-calcShrinkParameters(fitln,coef,smode,exclude)
    tauPos = shrinkPos$tau
    vpost = shrinkPos$v.post
    fitln[,"s2"] = vpost

    posRidge = sapply(seq(nr),function(i){
        k = which(posIndices[i,])
        y = log(mat[i,k])
        x = positiveMod[k,]
        l = vpost[i]/(nrow(x)*tauPos)
        if(i %in% exclude) return(matrix(rep(NA,ncol(positiveMod))))
        ridge = glmnet(y=y,x=x,lambda=l,alpha=0)
        as.matrix(coefficients(ridge)[colnames(positiveMod),])
    })
    posFittedCoefficients = t(posRidge)
    rownames(posFittedCoefficients) = rownames(mat)
    fitln[rownames(posFittedCoefficients),1:ncol(positiveMod)] = posFittedCoefficients
  }
  # calc zero component
  fitzero=calcZeroComponent(mat,zeroMod,posIndices)

  sdensity = density(fitzero[,"s2"],na.rm=TRUE)
  smode = sdensity$x[which.max(sdensity$y)]
  if(length(exclude)>0) fitzero[exclude,"s2"] = smode
  
  # shrink zero
  if(szero==TRUE){
    shrinkZero<-calcShrinkParameters(fitzero,coef,smode,exclude)
    tauZero = shrinkZero$tau
    vpostZero = shrinkZero$v.post
    fitzero[,"s2"] = vpostZero
    
    zeroRidge = sapply(1:nr,function(i){
      y = posIndices[i,]
      l = 1/(nc*tauZero)
      if(i %in% c(zeroExclude,exclude)) return(matrix(rep(NA,ncol(zeroMod))))
      ridge = glmnet(y=y,x=zeroMod,lambda=l,family="binomial",alpha=0,
        penalty.factor = c(rep(1,(ncol(zeroMod)-1)),0))
      as.matrix(coefficients(ridge))[colnames(zeroMod),]
    })
    zeroFittedCoefficients = t(zeroRidge)
    rownames(zeroFittedCoefficients) = rownames(mat)
    fitzero[rownames(zeroFittedCoefficients),1:ncol(zeroMod)] = zeroFittedCoefficients
  }

  # calc se
  se = calcStandardError(zeroMod,fitln,fitzero,coef=coef,exclude=union(exclude,zeroExclude))
  se[zeroExclude] = sqrt(fitln[zeroExclude,"s2"])

  # calc adjFactor
  adjFactor = calcZeroAdjustment(fitln,fitzero,zeroMod,coef,exclude=exclude)
  adjFactor[zeroExclude] = 0

  # calc logFC
  logFC <- fitln[,coef] + adjFactor

  list(logFC=logFC,adjFactor=adjFactor,se=se,
    fitln=fitln,fitzero=fitzero,zeroRidge=zeroRidge,posRidge=posRidge,
    tauPos=tauPos,tauZero=tauZero,exclude=exclude,zeroExclude=zeroExclude)
}
#' Positive component
#' 
#' Fit the positive (log-normal) component
#' 
#' @param mat A matrix of normalized counts
#' @param mod A model matrix
#' @param weights Weight matrix for samples and counts
#' @seealso \code{\link{fitZeroLogNormal}} \code{\link{fitFeatureModel}}
calcPosComponent<-function(mat,mod,weights){
  fitln <- lmFit(log(mat),mod,weights=weights)
  b = coefficients(fitln)
  df = fitln$df
  res = residuals(fitln,log(mat))
  s2 = sapply(seq(nrow(res)),function(i){
      sum(res[i,which(weights[i,])]^2,na.rm=TRUE)/df[i]
    })
  fitln<-data.frame(b=b,s2=s2,df=df)
  rownames(fitln) = rownames(mat)
  fitln
}
#' Zero component
#' 
#' Fit the zero (logisitic) component
#' 
#' @param mat A matrix of normalized counts
#' @param mod A model matrix
#' @param weights Weight matrix for samples and counts
#' @seealso \code{\link{fitZeroLogNormal}} \code{\link{fitFeatureModel}}
calcZeroComponent<-function(mat,mod,weights){
  fitzero <- sapply(seq(nrow(mat)), function(i) {
    fit <- glm.fit(mod, weights[i,], family=binomial())
    cf = coefficients(fit)
    df = fit$df.residual
    mc = exp(mod %*% cf)
    s2 = sum((weights[i, ] - t(mc/(1 + mc)))^2)/df    
    # s2 = sum(residuals(fit)^2)/df
    c(beta= cf, s2 = s2, df = df)
  })
  fitzero <- data.frame(t(fitzero))
  rownames(fitzero) = rownames(mat)
  fitzero
}
#' Calculate shrinkage parameters
#' 
#' Calculate the shrunken variances and variance of parameters of interest across features.
#' 
#' @param fit A matrix of fits as outputted by calcZeroComponent or calcPosComponent
#' @param coef Coefficient of interest
#' @param mins2 minimum variance estimate
#' @param exclude Vector of features to exclude when shrinking
#' @seealso \code{\link{fitZeroLogNormal}} \code{\link{fitFeatureModel}}
calcShrinkParameters<-function(fit,coef,mins2,exclude=NULL){

  if(is.null(exclude)){
    shrunkVar <- limma::squeezeVar(fit[,"s2"], fit[,"df"])
    v.post = shrunkVar$var.post
    tau <-var(fit[,coef],na.rm=TRUE)
  } else {
    v.post = rep(mins2,nrow(fit))
    shrunkVar <- limma::squeezeVar(fit[-exclude,"s2"], fit[-exclude,"df"])
    v.post[-exclude] <- shrunkVar$var.post
    tau <- var(fit[-exclude,coef],na.rm=TRUE)
  }
  list(tau=tau,v.post=v.post)
}
#' Calculate the zero-inflated component's adjustment factor
#' 
#' Calculate the log ratio of average marginal probabilities for each sample
#' having a positive count. This becomes the adjustment factor for the log 
#' fold change. 
#' 
#' @param fitln A matrix with parameters from the log-normal fit
#' @param fitzero A matrix with parameters from the logistic fit
#' @param mod The zero component model matrix
#' @param coef Coefficient of interest
#' @param exclude List of features to exclude
#' @seealso \code{\link{fitZeroLogNormal}} \code{\link{fitFeatureModel}}
calcZeroAdjustment<-function(fitln,fitzero,mod,coef,exclude=NULL){
  b = fitln[,1:(ncol(mod)-1)]
  beta = fitzero[,1:ncol(mod)]
  # calculate for zero adjust factor
  mod1 <- mod
  mod1[,coef] <- 1
  theta1 <- mod1 %*% t(beta)
  p1 <- exp(theta1) / (1+exp(theta1))
  p1 <- t(p1)
  if(ncol(b)>2) p1 = p1*exp(t(mod[,3:(ncol(mod)-1)]%*%t(b[,3:ncol(b)])))
  mean_p1 <- rowMeans(p1)

  mod0 <- mod
  mod0[,coef] <- 0
  theta0 <- mod0 %*% t(beta)
  p0 <- exp(theta0) / (1+exp(theta0))
  p0 <- t(p0)
  if(ncol(b)>2) p0 = p0*exp(t(mod[,3:(ncol(mod)-1)]%*%t(b[,3:ncol(b)])))
  mean_p0 <- rowMeans(p0)

  adjFactor <- log(mean_p1/mean_p0)
  if(!is.null(exclude)) adjFactor[exclude] = NA
  adjFactor
}

#' Calculate the zero-inflated log-normal statistic's standard error
#' 
#' Calculat the se for the model. Code modified from
#' "Adjusting for covariates in zero-inflated gamma and 
#' zero-inflated log-normal models for semicontinuous data", ED Mills 
#' 
#' @param mod The zero component model matrix
#' @param fitln A matrix with parameters from the log-normal fit
#' @param fitzero A matrix with parameters from the logistic fit
#' @param coef Coefficient of interest
#' @param exclude List of features to exclude
#' @seealso \code{\link{fitZeroLogNormal}} \code{\link{fitFeatureModel}}
calcStandardError<-function(mod,fitln,fitzero,coef=2,exclude=NULL){
  mod0 = mod1 = mod
  mod1[,coef] <- 1
  mod0[,coef] <- 0
  ve = rep(NA,nrow(fitln))
  features = seq(nrow(fitln))
  if(length(exclude)>0) features = features[-exclude]

# a) need to speed up
# b) need to include more covariates

  fullvar = sapply(features,function(i){
      beta = fitzero[i,1:ncol(mod)]
      b = fitln[i,1:(ncol(mod)-1)]
      s = as.numeric(fitln[i,"s2"])

      mu0 = as.vector(exp(mod0[,-ncol(mod)]%*%t(b) + .5*s))
      mu1 = as.vector(exp(mod1[,-ncol(mod)]%*%t(b) + .5*s))

      # calculate for zero adjust factor
      theta <- mod %*% t(beta)
      theta1 <- mod1 %*% t(beta)
      theta0 <- mod0 %*% t(beta)
      p  <- t(exp(theta) / (1+exp(theta)))
      p1 <- t(exp(theta1) / (1+exp(theta1)))
      p0 <- t(exp(theta0) / (1+exp(theta0)))

      checkInverse <- function(m){
        inherits(try(qr.solve(m),silent=T), "matrix")
      }
      
      Dp2 <- diag(length(p))*as.vector(p*(1-p))
      infz = t(mod)%*%Dp2%*%mod
      Dp <- diag(length(p))*as.vector(p)
      infln = t(mod[,-ncol(mod)])%*%Dp%*%mod[,-ncol(mod)]
      
      if(checkInverse(infz)) {
        invinf_z <-qr.solve(infz)
      } else {
        return(NA)
      }
      if(checkInverse(infln)) {
        invinf_ln<-as.numeric(s)*qr.solve(infln)
      } else {
        return(NA)
      }
      invInfFull = as.matrix( bdiag(invinf_z,invinf_ln, (2*s^2/sum(p))) )

      logRatioBeta0<- (mean(p1*(1-p1)*mu0)/mean(p1*mu0)) - (mean(p0*(1-p0)*mu0)/mean(p0*mu0))
      logRatioBeta1<-mean(p1*(1-p1)*mu0)/mean(p1*mu0)
      logRatioBeta2<- (mean(mod[,3]*p1*(1-p1)*mu0)/mean(p1*mu0)) - (mean(mod[,3]*p0*(1-p0)*mu0)/mean(p0*mu0))
      # logRatioB2<- (mean(mod[,3]*t(p1)*exp(mod0%*%t(b)))/mean(t(p1)*exp(mod0%*%t(b))))-
      #    (mean(mod[,3]*t(p0)*exp(mod0%*%t(b)))/mean(t(p0)*exp(mod0%*%t(b))))
      # logRatioFull = t(c(logRatioBeta0,logRatioBeta1,logRatioBeta2,0,1,logRatioB2,0))
      logRatioFull = t(c(logRatioBeta0,logRatioBeta1,logRatioBeta2,0,1,0))
      logRatioVar = logRatioFull%*%invInfFull%*%t(logRatioFull)
      logRatioVar
    })
  if(!is.null(exclude)){
    if(length(features)>0){
      ve[features] = fullvar
    }
  } else {
    ve = fullvar
  }
  sqrt(ve)
}
