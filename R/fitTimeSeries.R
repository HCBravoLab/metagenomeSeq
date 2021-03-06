#' Trapezoidal Integration
#' 
#' Compute the area of a function with values 'y' at the points 'x'.
#' Function comes from the pracma package.
#' 
#' @param x x-coordinates of points on the x-axis
#' @param y y-coordinates of function values
#' @return Approximated integral of the function from 'min(x)' to 'max(x)'. 
#'  Or a matrix of the same size as 'y'.
#' @rdname trapz
#' @export
#' @examples
#' 
#' # Calculate the area under the sine curve from 0 to pi:
#'  n <- 101
#'  x <- seq(0, pi, len = n)
#'  y <- sin(x)
#'  trapz(x, y)          #=> 1.999835504
#' 
#' # Use a correction term at the boundary: -h^2/12*(f'(b)-f'(a))
#'  h  <- x[2] - x[1]
#'  ca <- (y[2]-y[1]) / h
#'  cb <- (y[n]-y[n-1]) / h
#'  trapz(x, y) - h^2/12 * (cb - ca)  #=> 1.999999969
#'
trapz <- function(x,y){
    if (missing(y)) {
        if (length(x) == 0) 
            return(0)
        y <- x
        x <- 1:length(x)
    }
    if (length(x) == 0) 
        return(0)
    if (!(is.numeric(x) || is.complex(x)) || !(is.numeric(y) || 
        is.complex(y))) 
        stop("Arguments 'x' and 'y' must be real or complex.")
    m <- length(x)
    xp <- c(x, x[m:1])
    yp <- c(numeric(m), y[m:1])
    n <- 2 * m
    p1 <- sum(xp[1:(n - 1)] * yp[2:n]) + xp[n] * yp[1]
    p2 <- sum(xp[2:n] * yp[1:(n - 1)]) + xp[1] * yp[n]
    return(0.5 * (p1 - p2))
}

#' smoothing-splines anova fit
#' 
#' Sets up a data-frame with the feature abundance, 
#' class information, time points, sample ids and returns
#' the fitted values for the fitted model.
#' 
#' @param formula Formula for ssanova. Of the form: abundance ~ ... where ... includes any pData slot value.
#' @param abundance Numeric vector of abundances.
#' @param class Class membership (factor of group membership).
#' @param time Time point vector of relative times (same length as abundance).
#' @param id Sample / patient id.
#' @param include Parameters to include in prediction.
#' @param pd Extra variable.
#' @param ... Extra parameters for ssanova function (see ?ssanova).
#' @return \itemize{A list containing:
#' \item     data        : Inputed data
#' \item     fit         : The interpolated / fitted values for timePoints
#' \item     se          : The standard error for CI intervals
#' \item     timePoints  : The time points interpolated over
#' }
#' @seealso \code{\link{cumNorm}} \code{\link{fitTimeSeries}} \code{\link{ssPermAnalysis}} \code{\link{ssPerm}} \code{\link{ssIntervalCandidate}}
#' @rdname ssFit
#' @export
#' @examples
#'
#' # Not run
#'
ssFit <- function(formula,abundance,class,time,id,include=c("class", "time:class"),pd,...) {
    df = data.frame(abundance = abundance, class = factor(class),
       time=time,id = factor(id),pd)
    
    # The smoothing splines anova model
    if(missing(formula)){
        mod = gss::ssanova(abundance ~ time * class, data=df,...)
    } else{
        mod = gss::ssanova(formula,data=df,...)
    }

    fullTime = seq(min(df$time), max(df$time), by=1)
    values = data.frame(time=fullTime, class=factor(levels(df[,"class"]))[2])
    fit = predict(mod, values, include=include, se=TRUE)
    
    res = list(data=df, fit=fit$fit, se=fit$se, timePoints=fullTime)
    return(res)
}

#' class permutations for smoothing-spline time series analysis
#' 
#' Creates a list of permuted class memberships for the time series permuation tests.
#' 
#' @param df Data frame containing class membership and sample/patient id label.
#' @param B Number of permutations.
#' @return A list of permutted class memberships
#' @seealso \code{\link{cumNorm}} \code{\link{fitTimeSeries}} \code{\link{ssFit}} \code{\link{ssPermAnalysis}} \code{\link{ssIntervalCandidate}}
#' @rdname ssPerm
#' @examples
#'
#' # Not run
#'
ssPerm <- function(df,B) {
    dat = data.frame(class=df$class, id=df$id)
    # id  = table(dat$id)
    id = table(interaction(dat$class,dat$id))
    id = id[id>0]
    classes = unique(dat)[,"class"]
    permList = lapply(1:B,function(i){
        rep(sample(classes, replace=FALSE),id)
    })
    return(permList) 
}

#' smoothing-splines anova fits for each permutation
#' 
#' Calculates the fit for each permutation and estimates 
#' the area under the null (permutted) model for interesting time 
#' intervals of differential abundance.
#' 
#' @param data Data used in estimation.
#' @param formula Formula for ssanova. Of the form: abundance ~ ... where ... includes any pData slot value.
#' @param permList A list of permutted class memberships
#' @param intTimes Interesting time intervals.
#' @param timePoints Time points to interpolate over.
#' @param include Parameters to include in prediction.
#' @param ... Options for ssanova
#' @return A matrix of permutted area estimates for time intervals of interest.
#' @seealso \code{\link{cumNorm}} \code{\link{fitTimeSeries}} \code{\link{ssFit}} \code{\link{ssPerm}} \code{\link{ssIntervalCandidate}}
#' @rdname ssPermAnalysis
#' @export
#' @examples
#'
#' # Not run
#'
ssPermAnalysis <- function(data,formula,permList,intTimes,timePoints,include=c("class", "time:class"),...){
    resPerm=matrix(NA, length(permList), nrow(intTimes))
    permData=data
    case = data.frame(time=timePoints, class=factor(levels(data$class)[2]))
    for (j in 1:length(permList)){
        
        permData$class = permList[[j]]
        # The smoothing splines anova model
        if(!missing(formula)){
            permModel = gss::ssanova(formula, data=permData,...)
        } else{
            permModel = gss::ssanova(abundance ~ time * class,data=permData,...)
        }

        permFit = cbind(timePoints, (2*predict(permModel,case,include=include, se=TRUE)$fit))
            for (i in 1:nrow(intTimes)){
                permArea=permFit[which(permFit[,1]==intTimes[i,1]) : which(permFit[,1]==intTimes[i, 2]), ]
                resPerm[j, i]=metagenomeSeq::trapz(x=permArea[,1], y=permArea[,2])
            }
        if(j%%100==0) show(j)
    }
    return(resPerm)
}

#' calculate interesting time intervals
#' 
#' Calculates time intervals of interest using SS-Anova fitted confidence intervals.
#' 
#' @param fit SS-Anova fits.
#' @param standardError SS-Anova se estimates.
#' @param timePoints Time points interpolated over.
#' @param positive Positive region or negative region (difference in abundance is positive/negative).
#' @param C Value for which difference function has to be larger or smaller than (default 0).
#' @return Matrix of time point intervals of interest
#' @seealso \code{\link{cumNorm}} \code{\link{fitTimeSeries}} \code{\link{ssFit}} \code{\link{ssPerm}} \code{\link{ssPermAnalysis}}
#' @rdname ssIntervalCandidate
#' @export
#' @examples
#'
#' # Not run
#'
ssIntervalCandidate <- function(fit, standardError, timePoints, positive=TRUE,C=0){
    lowerCI = (2*fit - (1.96*2*standardError))
    upperCI = (2*fit + (1.96*2*standardError))
    if (positive){
        abundanceDifference = which( lowerCI>=0 & abs(lowerCI)>=C )
    }else{
        abundanceDifference = which( upperCI<=0 & abs(upperCI)>=C )
    }
    if (length(abundanceDifference)>0){
        intIndex=which(diff(abundanceDifference)!=1)
        intTime=matrix(NA, (length(intIndex)+1), 4)
        if (length(intIndex)==0){
            intTime[1,1]=timePoints[abundanceDifference[1]]
            intTime[1,2]=timePoints[tail(abundanceDifference, n=1)]
        }else{
            i=1
            while(length(intTime)!=0 & length(intIndex)!=0){
                intTime[i,1]=timePoints[abundanceDifference[1]]
                intTime[i,2]=timePoints[abundanceDifference[intIndex[1]]]
                abundanceDifference=abundanceDifference[-c(1:intIndex[1])]
                intIndex=intIndex[-1]
                i=i+1
            }
        intTime[i,1] = timePoints[abundanceDifference[1]]
        intTime[i,2] = timePoints[tail(abundanceDifference, n=1)]
        }
    }else{
        intTime=NULL   
    }   
    return(intTime)    
}

#' Discover differentially abundant time intervals using SS-Anova
#' 
#' Calculate time intervals of interest using SS-Anova fitted models.
#' Fitting is performed uses Smoothing Spline ANOVA (SS-Anova) to find interesting intervals of time. 
#' Given observations at different time points for two groups, fitSSTimeSeries 
#' calculates a  function that models the difference in abundance between two 
#' groups across all time. Using permutations we estimate a null distribution 
#' of areas for the time intervals of interest and report significant intervals of time.
#' Use of the function for analyses should cite:
#' "Finding regions of interest in high throughput genomics data using smoothing splines"
#' Talukder H, Paulson JN, Bravo HC. (In preparation)
#' 
#' @param obj metagenomeSeq MRexperiment-class object.
#' @param formula Formula for ssanova. Of the form: abundance ~ ... where ... includes any pData slot value.
#' @param feature Name or row of feature of interest.
#' @param class Name of column in phenoData of MRexperiment-class object for class memberhip.
#' @param time Name of column in phenoData of MRexperiment-class object for relative time.
#' @param id Name of column in phenoData of MRexperiment-class object for sample id.
#' @param lvl Vector or name of column in featureData of MRexperiment-class object for aggregating counts (if not OTU level).
#' @param include Parameters to include in prediction.
#' @param C Value for which difference function has to be larger or smaller than (default 0).
#' @param B Number of permutations to perform
#' @param norm When aggregating counts to normalize or not.
#' @param log Log2 transform.
#' @param sl Scaling value.
#' @param featureOrder Hierarchy of levels in taxonomy as fData colnames
#' @param ... Options for ssanova
#' @return List of matrix of time point intervals of interest, Difference in abundance area and p-value, fit, area permutations, and call.
#' @return A list of objects including:
#' \itemize{
#'  \item{timeIntervals - Matrix of time point intervals of interest, area of differential abundance, and pvalue.}
#'  \item{data  - Data frame of abundance, class indicator, time, and id input.}
#'  \item{fit - Data frame of fitted values of the difference in abundance, standard error estimates and timepoints interpolated over.}
#'  \item{perm - Differential abundance area estimates for each permutation.}
#'  \item{call - Function call.}
#' }
#' @rdname fitSSTimeSeries
#' @seealso \code{\link{cumNorm}} \code{\link{ssFit}} \code{\link{ssIntervalCandidate}} \code{\link{ssPerm}} \code{\link{ssPermAnalysis}} \code{\link{plotTimeSeries}}
#' @export
#' @examples
#'
#' data(mouseData)
#' res = fitSSTimeSeries(obj=mouseData,feature="Actinobacteria",
#'    class="status",id="mouseID",time="relativeTime",lvl='class',B=2)
#'
fitSSTimeSeries <- function(obj,formula,feature,class,time,id,lvl=NULL,include=c("class", "time:class"),C=0,B=1000,norm=TRUE,log=TRUE,sl=1000,featureOrder=NULL,...) {
    
    if(!is.null(lvl)){
        aggData = aggregateByTaxonomy(obj,lvl,norm=norm,sl=sl, featureOrder=featureOrder)
        abundance = MRcounts(aggData,norm=FALSE,log=log,sl=1)[feature,]
    } else { 
        abundance = MRcounts(obj,norm=norm,log=log,sl=sl)[feature,]
    }
    class = pData(obj)[,class]
    time  = pData(obj)[,time]
    id    = pData(obj)[,id]
    if(any(sapply(list(id,time,class),length)==0)){
        stop("provide class, time, and id names")
    }

    if(!missing(formula)){
        prep=ssFit(formula=formula,abundance=abundance,class=class,
            time=time,id=id,include=include,pd=pData(obj),...)
    } else {
        prep=ssFit(abundance=abundance,class=class,time=time,id=id,
            include=include,pd=pData(obj),...)
    }
    indexPos = ssIntervalCandidate(fit=prep$fit, standardError=prep$se, 
        timePoints=prep$timePoints, positive=TRUE,C=C)
    indexNeg = ssIntervalCandidate(fit=prep$fit, standardError=prep$se, 
        timePoints=prep$timePoints, positive=FALSE,C=C)
    indexAll = rbind(indexPos, indexNeg)

    if(sum(indexAll[,1]==indexAll[,2])>0){
        indexAll=indexAll[-which(indexAll[,1]==indexAll[,2]),]
    }

    fit = 2*prep$fit
    se  = 2*prep$se
    timePoints = prep$timePoints
    fits = data.frame(fit = fit, se = se, timePoints = timePoints)
    
    if(!is.null(indexAll)){
      if(length(indexAll)>0){
        indexAll=matrix(indexAll,ncol=4)
        colnames(indexAll)=c("Interval start", "Interval end", "Area", "p.value")
        predArea    = cbind(prep$timePoints, (2*prep$fit))
        permList    = ssPerm(prep$data,B=B)
        if(!missing(formula)){
            permResult  = ssPermAnalysis(data=prep$data,formula=formula,permList=permList,
                intTimes=indexAll,timePoints=prep$timePoints,include=include,...)
        } else {
            permResult  = ssPermAnalysis(data=prep$data,permList=permList,
                intTimes=indexAll,timePoints=prep$timePoints,include=include,...)
        }
        
        for (i in 1:nrow(indexAll)){
            origArea=predArea[which(predArea[,1]==indexAll[i,1]):which(predArea[,1]==indexAll[i, 2]), ]
            actArea=trapz(x=origArea[,1], y=origArea[,2])
            indexAll[i,3] = actArea
            if(actArea>0){
                indexAll[i,4] = 1 - (length(which(actArea>permResult[,i]))+1)/(B+1)
            }else{
                indexAll[i,4] = (length(which(actArea>permResult[,i]))+1)/(B+1)
            }
        if(indexAll[i,4]==0){ 
        indexAll[i,4] = 1/(B+1)
        }
        }

        res = list(timeIntervals=indexAll,data=prep$data,fit=fits,perm=permResult)
        return(res)
      }
    }else{
        indexAll = "No statistically significant time intervals detected"
        res = list(timeIntervals=indexAll,data=prep$data,fit=fits,perm=NULL)
        return(res)
    }
}

#' Discover differentially abundant time intervals
#' 
#' Calculate time intervals of significant differential abundance.
#' Currently only one method is implemented (ssanova). fitSSTimeSeries is called with method="ssanova".
#' 
#' @param obj metagenomeSeq MRexperiment-class object.
#' @param formula Formula for ssanova. Of the form: abundance ~ ... where ... includes any pData slot value.
#' @param feature Name or row of feature of interest.
#' @param class Name of column in phenoData of MRexperiment-class object for class memberhip.
#' @param time Name of column in phenoData of MRexperiment-class object for relative time.
#' @param id Name of column in phenoData of MRexperiment-class object for sample id.
#' @param method Method to estimate time intervals of differentially abundant bacteria (only ssanova method implemented currently).
#' @param lvl Vector or name of column in featureData of MRexperiment-class object for aggregating counts (if not OTU level).
#' @param include Parameters to include in prediction.
#' @param C Value for which difference function has to be larger or smaller than (default 0).
#' @param B Number of permutations to perform.
#' @param norm When aggregating counts to normalize or not.
#' @param log Log2 transform.
#' @param sl Scaling value.
#' @param featureOrder Hierarchy of levels in taxonomy as fData colnames
#' @param ... Options for ssanova
#' @return List of matrix of time point intervals of interest, Difference in abundance area and p-value, fit, area permutations, and call.
#' @return A list of objects including:
#' \itemize{
#'  \item{timeIntervals - Matrix of time point intervals of interest, area of differential abundance, and pvalue.}
#'  \item{data  - Data frame of abundance, class indicator, time, and id input.}
#'  \item{fit - Data frame of fitted values of the difference in abundance, standard error estimates and timepoints interpolated over.}
#'  \item{perm - Differential abundance area estimates for each permutation.}
#'  \item{call - Function call.}
#' }
#' @rdname fitTimeSeries
#' @seealso \code{\link{cumNorm}} \code{\link{fitSSTimeSeries}} \code{\link{plotTimeSeries}}
#' @export
#' @examples
#'
#' data(mouseData)
#' res = fitTimeSeries(obj=mouseData,feature="Actinobacteria",
#'    class="status",id="mouseID",time="relativeTime",lvl='class',B=2)
#'
fitTimeSeries <- function(obj,formula,feature,class,time,id,method=c("ssanova"),
                        lvl=NULL,include=c("class", "time:class"),C=0,B=1000,
                        norm=TRUE,log=TRUE,sl=1000,featureOrder=NULL,...) {
    if(method=="ssanova"){
        if(requireNamespace("gss")){
            if(missing(formula)){
                res = fitSSTimeSeries(obj=obj,feature=feature,class=class,time=time,id=id,
                        lvl=lvl,C=C,B=B,norm=norm,log=log,sl=sl,include=include,featureOrder=featureOrder,...)
            } else {
                res = fitSSTimeSeries(obj=obj,formula=formula,feature=feature,class=class,
                        time=time,id=id,lvl=lvl,C=C,B=B,norm=norm,log=log,sl=sl,
                        include=include,featureOrder=featureOrder,...)
            }
        }
    }
    res = c(res,call=match.call())
    return(res)
}

#' Plot difference function for particular bacteria
#' 
#' Plot the difference in abundance for significant features.
#' 
#' @param res Output of fitTimeSeries function
#' @param C Value for which difference function has to be larger or smaller than (default 0).
#' @param xlab X-label.
#' @param ylab Y-label.
#' @param main Main label.
#' @param ... Extra plotting arguments.
#' @return Plot of difference in abundance for significant features.
#' @rdname plotTimeSeries
#' @seealso \code{\link{fitTimeSeries}}
#' @export
#' @examples
#'
#' data(mouseData)
#' res = fitTimeSeries(obj=mouseData,feature="Actinobacteria",
#'    class="status",id="mouseID",time="relativeTime",lvl='class',B=10)
#' plotTimeSeries(res)
#'
plotTimeSeries<-function(res,C=0,xlab="Time",ylab="Difference in abundance",main="SS difference function prediction",...){
    fit = res$fit$fit
    se  = res$fit$se
    timePoints = res$fit$timePoints
    confInt95 = 1.96
    sigDiff = res$timeIntervals

    minValue=min(fit-(confInt95*se))-.5
    maxValue=max(fit+(confInt95*se))+.5

    plot(x=timePoints, y=fit, ylim=c(minValue, maxValue), xlab=xlab, ylab=ylab, main=main, ...)

    for (i in 1:nrow(sigDiff)){
        begin=sigDiff[i,1]
        end=sigDiff[i,2]
        indBegin=which(timePoints==begin)
        indEnd=which(timePoints==end)
        x=timePoints[indBegin:indEnd]
        y=fit[indBegin:indEnd]
        xx=c(x, rev(x))
        yy=c(y, rep(0, length(y)))
        polygon(x=xx, yy, col="grey")
    }
    lines(x=timePoints, y=fit, pch="")
    lines(x=timePoints, y=fit+(confInt95*se), pch="", lty=2)
    lines(x=timePoints, y=fit-(confInt95*se), pch="", lty=2)
    abline(h=C)
}

#' Plot abundances by class
#' 
#' Plot the abundance of values for each class using 
#' a spline approach on the estimated full model.
#' 
#' @param res Output of fitTimeSeries function
#' @param formula Formula for ssanova. Of the form: abundance ~ ... where ... includes any pData slot value.
#' @param xlab X-label.
#' @param ylab Y-label.
#' @param color0 Color of samples from first group.
#' @param color1 Color of samples from second group.
#' @param include Parameters to include in prediction.
#' @param ... Extra plotting arguments.
#' @return Plot for abundances of each class using a spline approach on estimated null model.
#' @rdname plotClassTimeSeries
#' @seealso \code{\link{fitTimeSeries}}
#' @export
#' @examples
#'
#' data(mouseData)
#' res = fitTimeSeries(obj=mouseData,feature="Actinobacteria",
#'    class="status",id="mouseID",time="relativeTime",lvl='class',B=10)
#' plotClassTimeSeries(res,pch=21,bg=res$data$class,ylim=c(0,8))
#'
plotClassTimeSeries<-function(res,formula,xlab="Time",ylab="Abundance",color0="black",
                            color1="red",include=c("1","class", "time:class"),...){
    dat = res$data
    if(missing(formula)){
        mod = gss::ssanova(abundance ~ time * class, data=dat)
    } else{
        mod = gss::ssanova(formula,data=dat)
    }
    
    timePoints = seq(min(dat$time),max(dat$time),by=1)
    group0 = data.frame(time=timePoints,class=factor(levels(dat$class)[1]))
    group1 = data.frame(time=timePoints,class=factor(levels(dat$class)[2]))

    pred0  = predict(mod, newdata=group0,include=include, se=TRUE)
    pred1  = predict(mod, newdata=group1,include=include, se=TRUE)
    
    plot(x=dat$time,y=dat$abundance,xlab=xlab,ylab=ylab,...)
    lines(x=group0$time,y=pred0$fit,col=color0)
    lines(x=group0$time,y=pred0$fit+(1.96*pred0$se),lty=2,col=color0)
    lines(x=group0$time,y=pred0$fit-(1.96*pred0$se),lty=2,col=color0)

    lines(x=group1$time,y=pred1$fit,col=color1)
    lines(x=group1$time,y=pred1$fit+(1.96*pred1$se),lty=2,col=color1)
    lines(x=group1$time,y=pred1$fit-(1.96*pred1$se),lty=2,col=color1)
}

#' Discover differentially abundant time intervals for all bacteria
#' 
#' Calculate time intervals of significant differential abundance over all
#' bacteria of a particularly specified level (lvl). If not lvl is specified,
#' all OTUs are analyzed. Warning, function can take a while
#' 
#' @param obj metagenomeSeq MRexperiment-class object.
#' @param lvl Vector or name of column in featureData of MRexperiment-class object for aggregating counts (if not OTU level).
#' @param B Number of permutations to perform.
#' @param featureOrder Hierarchy of levels in taxonomy as fData colnames
#' @param ... Options for \code{\link{fitTimeSeries}}, except feature.
#' @return List of lists of matrices of time point intervals of interest, Difference in abundance area and p-value, fit, area permutations.
#' @return A list of lists for which each includes:
#' \itemize{
#'  \item{timeIntervals - Matrix of time point intervals of interest, area of differential abundance, and pvalue.}
#'  \item{data  - Data frame of abundance, class indicator, time, and id input.}
#'  \item{fit - Data frame of fitted values of the difference in abundance, standard error estimates and timepoints interpolated over.}
#'  \item{perm - Differential abundance area estimates for each permutation.}
#'  \item{call - Function call.}
#' }
#' @rdname fitMultipleTimeSeries
#' @seealso \code{\link{cumNorm}} \code{\link{fitSSTimeSeries}} \code{\link{fitTimeSeries}}
#' @export
#' @examples
#'
#' data(mouseData)
#' res = fitMultipleTimeSeries(obj=mouseData,lvl='phylum',class="status",
#'           id="mouseID",time="relativeTime",B=1)
#'
fitMultipleTimeSeries <- function(obj,lvl=NULL,B=1,featureOrder=NULL,...) {
    if(is.null(lvl)){
        bacteria = seq(nrow(obj))
    } else {
        if(is.factor(fData(obj)[,lvl])){
            fData(obj)[,lvl] = as.character(fData(obj)[,lvl])
        }
        bacteria = unique(fData(obj)[,lvl])
    }
    fits = lapply(bacteria,function(bact){
        try(fitTimeSeries(obj,lvl=lvl,feature=bact,B=B,featureOrder=featureOrder,...))
    })
    names(fits) = bacteria
    fits = c(fits,call=match.call())
    return(fits)
}

#' With a list of fitTimeSeries results, generate
#' an MRexperiment that can be plotted with metavizr
#' 
#' @param obj Output of fitMultipleTimeSeries
#' @param sampleNames Sample names for plot
#' @param sampleDescription Description of samples for plot axis label
#' @param taxonomyLevels Feature names for plot
#' @param taxonomyHierarchyRoot Root of feature hierarchy for MRexperiment
#' @param taxonomyDescription Description of features for plot axis label
#' @param featuresOfInterest The features to select from the fitMultipleTimeSeries output
#' @param featureDataOfInterest featureData for the resulting MRexperiment
#' @return MRexperiment that contains fitTimeSeries data, featureData, and phenoData
#' @rdname ts2MRexperiment
#' @seealso \code{\link{fitTimeSeries}} \code{\link{fitMultipleTimeSeries}}
#' @export
#' @examples
#'
#' data(mouseData)
#' res = fitMultipleTimeSeries(obj=mouseData,lvl='phylum',class="status",
#'           id="mouseID",time="relativeTime",B=1)
#' obj = ts2MRexperiment(res)
#' obj
#'
ts2MRexperiment<-function(obj,sampleNames=NULL,
                          sampleDescription="timepoints",
                          taxonomyLevels=NULL,
                          taxonomyHierarchyRoot="bacteria",
                          taxonomyDescription="taxonomy",
                          featuresOfInterest = NULL,
                          featureDataOfInterest=NULL){
  if(is.null(obj)){
    stop("Matrix cannot be null")
  }
  if(is.null(sampleNames)){
    numSamples <- dim(obj[[1]]$fit)[1]
    sampleNames <- paste("Timepoint", 1:numSamples, sep="_")
  }
  
  if(is.null(featuresOfInterest)){
    hasFit <- lapply(1:(length(obj)-1), function(i) which(!is.null(obj[[i]]$fit)))
    featuresOfInterest <- which(hasFit == 1)
    hasFit <- (hasFit == 1)
    hasFit <- !is.na(hasFit)
    temp <- 1:length(hasFit)
    temp[!hasFit] <- 0
    hasFit <- temp
  }
  
  if(is.null(taxonomyLevels)){
    numLevels <- 1:length(hasFit)
    taxonomyLevels <- names(obj)[1:length(hasFit)]
  }
  
  numSamples <- length(sampleNames)
  numLevels <- length(taxonomyLevels)
  numFeaturesOfInterest <- length(featuresOfInterest)
  
  rangeSamples <- 1:numSamples
  rangeFeaturesOfInterest <- 1:numFeaturesOfInterest
#  print(hasFit)
  
  results <- do.call(rbind, lapply(hasFit,function(i){ if (i != 0) t(obj[[i]]$fit)[1,] else rep(NA, numSamples) }))
  
  dfSamples <- data.frame(x=rangeSamples,row.names=sampleNames)
  metaDataSamples <-data.frame(labelDescription=sampleDescription)
  annotatedDFSamples <- AnnotatedDataFrame()
  pData(annotatedDFSamples) <- dfSamples
  varMetadata(annotatedDFSamples) <- metaDataSamples
  validObject(annotatedDFSamples)
  
  if(is.null(featureDataOfInterest)){
    dfFeatures <- data.frame(taxonomy1=rep(taxonomyHierarchyRoot, numLevels),taxonomy2=taxonomyLevels)
    metaDataFeatures <-data.frame(labelDescription=paste(taxonomyDescription, 1:2, sep=""))
    annotatedDFFeatures <- AnnotatedDataFrame()
    pData(annotatedDFFeatures) <- dfFeatures
    varMetadata(annotatedDFFeatures) <- metaDataFeatures
    validObject(annotatedDFFeatures)
  }
  else{
    annotatedDFFeatures <- featureDataOfInterest
  }
  
  fitTimeSeriesMRexp <- newMRexperiment(counts=results,
                                        phenoData=annotatedDFSamples,
                                        featureData=annotatedDFFeatures)
  return(fitTimeSeriesMRexp)
}
# load("~/Dropbox/Projects/metastats/package/git/metagenomeSeq/data/mouseData.rda")
# classMatrix = aggregateByTaxonomy(mouseData,lvl='class',norm=TRUE,out='MRexperiment')
# data(mouseData)
# fitTimeSeries(obj=mouseData,feature="Actinobacteria",class="status",id="mouseID",time="relativeTime",lvl='class',B=10)
