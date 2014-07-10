#' @name trapz
#' @title Trapezoidal Integration
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

#' @name ssFit
#' @title smoothing-splines anova fit
#' 
#' @details Sets up a data-frame with the feature abundance, 
#' class information, time points, sample ids and returns
#' the fitted values for the fitted model.
#' 
#' @param abundance Numeric vector of abundances.
#' @param class Class membership (factor of group membership).
#' @param time Time point vector of relative times (same length as abundance).
#' @param id Sample / patient id.
#' @param ... Extra parameters for ssanova function (see ?ssanova).
#' @return A list containing:
#'      data        : Inputed data
#'      fit         : The interpolated / fitted values for timePoints
#'      se          : The standard error for CI intervals
#'      timePoints  : The time points interpolated over
#' @seealso \code{\link{cumNorm}} \code{\link{fitTimeSeries}} \code{\link{ssPermAnalysis}} \code{\link{ssPerm}} \code{\link{ssIntervalCandidate}}
#' @rdname ssFit
#' @export
#' @examples
#'
#' # Not run
#'
ssFit <- function(abundance,class,time,id,...) {
    df = data.frame(abundance = abundance, class = factor(class),
        time=time,id = factor(id))

    # The smoothing splines anova model
    mod = ssanova(abundance ~ time * class, data=df,...)
    fullTime = seq(min(df$time), max(df$time), by=1)
    values = data.frame(time=fullTime, class=factor(levels(df[,"class"]))[1])
    fit = predict(mod, values, include=c("class", "time:class"), se=TRUE)
    
    res = list(data=df, fit=fit$fit, se=fit$se, timePoints=fullTime)
    return(res)
}

#' @name ssPerm
#' @title class permutations for smoothing-spline time series analysis
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
    id  = table(dat$id)
    classes = unique(dat)[,"class"]
    permList = lapply(1:B,function(i){
        rep(sample(classes, replace=FALSE),id)
    })
    return(permList) 
}

#' @name ssPermAnalysis
#' @title smoothing-splines anova fits for each permutation
#' 
#' @details Calculates the fit for each permutation and estimates 
#' the area under the null (permutted) model for interesting time 
#' intervals of differential abundance.
#' 
#' @param data Data used in estimation.
#' @param permList A list of permutted class memberships
#' @param intTimes Interesting time intervals.
#' @param timePoints Time points to interpolate over.
#' @param ... Options for ssanova
#' @return A matrix of permutted area estimates for time intervals of interest.
#' @seealso \code{\link{cumNorm}} \code{\link{fitTimeSeries}} \code{\link{ssFit}} \code{\link{ssPerm}} \code{\link{ssIntervalCandidate}}
#' @rdname ssPermAnalysis
#' @export
#' @examples
#'
#' # Not run
#'
ssPermAnalysis <- function(data, permList, intTimes, timePoints,...){
    resPerm=matrix(NA, length(permList), nrow(intTimes))
    permData=data
    for (j in 1:length(permList)){
        
        permData$class = permList[[j]]
        permModel      = ssanova(abundance ~ time * class, data=permData,...)
        permFit        = cbind(timePoints, (2*predict(permModel,data.frame(time=timePoints, class=factor(1)),#abs 
            include=c("class", "time:class"), se=TRUE)$fit))

            for (i in 1:nrow(intTimes)){
                permArea=permFit[which(permFit[,1]==intTimes[i,1]) : which(permFit[,1]==intTimes[i, 2]), ]
                resPerm[j, i]=metagenomeSeq::trapz(x=permArea[,1], y=permArea[,2])
            }
        if(j%%100==0) show(j)
    }
    return(resPerm)
}

#' @name ssIntervalCandidate
#' @title calculate interesting time intervals
#' 
#' Calculates time intervals of interest using SS-Anova fitted confidence intervals.
#' 
#' @param fit SS-Anova fits.
#' @param standardError SS-Anova se estimates.
#' @param timePoints Time points interpolated over.
#' @param positive Positive region or negative region (difference in abundance is positive/negative).
#' @param C Value for which difference function has to be larger than (default 0).
#' @return Matrix of time point intervals of interest
#' @seealso \code{\link{cumNorm}} \code{\link{fitTimeSeries}} \code{\link{ssFit}} \code{\link{ssPerm}} \code{\link{ssPermAnalysis}}
#' @rdname ssIntervalCandidate
#' @export
#' @examples
#'
#' # Not run
#'
ssIntervalCandidate <- function(fit, standardError, timePoints, positive=TRUE,C=0){
    if (positive){
        abundanceDifference = which( (2*fit - (1.96*2*standardError) )>C)
    }else{
        abundanceDifference = which( (2*fit + (1.96*2*standardError) )<C)
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

#' @name fitTimeSeries
#' @title Discover differentially abundant time intervals using SS-Anova
#' 
#' @details Calculate time intervals of interest using SS-Anova fitted models.
#' Fitting is performed uses Smoothing Spline ANOVA (SS-Anova) to find interesting intervals of time. 
#' Given observations at different time points for two groups, fitTimeSeries 
#' calculates a  function that models the difference in abundance between two 
#' groups across all time. Using permutations we estimate a null distribution 
#' of areas for the time intervals of interest and report significant intervals of time.
#' Use of the function for analyses should cite:
#' "Finding regions of interest in high throughput genomics data using smoothing splines"
#' Talukder H, Paulson JN, Bravo HC. (Submitted)
#' 
#' @param obj metagenomeSeq MRexperiment-class object.
#' @param feature Name or row of feature of interest.
#' @param class Name of column in phenoData of MRexperiment-class object for class memberhip.
#' @param time Name of column in phenoData of MRexperiment-class object for relative time.
#' @param id Name of column in phenoData of MRexperiment-class object for sample id.
#' @param lvl Vector or name of column in featureData of MRexperiment-class object for aggregating counts (if not OTU level).
#' @param C Value for which difference function has to be larger than.
#' @param B Number of permutations to perform
#' @param seed Random-number seed.
#' @param norm When aggregating counts to normalize or not.
#' @param sl Scaling value.
#' @param ... Options for ssanova
#' @return List of matrix of time point intervals of interest, Difference in abundance area and p-value, fit, area permutations, and call.
#' @rdname fitTimeSeries
#' @seealso \code{\link{cumNorm}} \code{\link{ssFit}} \code{\link{ssIntervalCandidate}} \code{\link{ssPerm}} \code{\link{ssPermAnalysis}} \code{\link{plotTimeSeries}}
#' @export
#' @examples
#'
#' data(mouseData)
#' res = fitTimeSeries(obj=mouseData,feature="Actinobacteria",
#'    class="status",id="mouseID",time="relativeTime",lvl='class',B=10)
#'
fitTimeSeries <- function(obj,feature,class,time,id,lvl=NULL,C=0,B=1000,seed=123,norm=TRUE,sl=1000,...) {
    if(!require(gss)){
        install.packages("gss",repos="http://cran.r-project.org")
        library(gss)
    }
    set.seed(seed)
    
    if(!is.null(lvl)){
        aggData = aggregateByTaxonomy(obj,lvl,norm=norm,sl=sl)
        abundance = MRcounts(aggData,norm=FALSE,log=TRUE,sl=1)[feature,]
    } else { 
        abundance = MRcounts(obj,norm=norm,log=TRUE,sl=sl)
    }
    class = pData(obj)[,class]
    time  = pData(obj)[,time]
    id    = pData(obj)[,id]

    prep=ssFit(abundance=abundance,class=class,time=time,id=id,...)
    indexPos = ssIntervalCandidate(fit=prep$fit, standardError=prep$se, 
        timePoints=prep$timePoints, positive=TRUE,C=C)
    indexNeg = ssIntervalCandidate(fit=prep$fit, standardError=prep$se, 
        timePoints=prep$timePoints, positive=FALSE,C=C)
    indexAll = rbind(indexPos, indexNeg)

    if(!is.null(indexAll)){
        colnames(indexAll)=c("Interval start", "Interval end", "Area", "p.value")
        predArea    = cbind(prep$timePoints, (2*prep$fit))
        permList    = ssPerm(prep$data,B=B)
        permResult  = ssPermAnalysis(data=prep$data, permList=permList,
            intTimes=indexAll, timePoints=prep$timePoints,...)
        
        for (i in 1:nrow(indexAll)){
            origArea=predArea[which(predArea[,1]==indexAll[i,1]):which(predArea[,1]==indexAll[i, 2]), ]
            actArea=trapz(x=origArea[,1], y=origArea[,2])
            indexAll[i,3] = actArea
            if(actArea>0){
                indexAll[i,4] = 1 - length(which(actArea>permResult[,i]))/B
            }else{
                indexAll[i,4] = length(which(actArea>permResult[,i]))/B
            }

        }
        fit = 2*prep$fit
        se  = 2*prep$se
        timePoints = prep$timePoints
        fits = data.frame(fit = fit, se = se, timePoints = timePoints)

        res = list(timeIntervals=indexAll,data=prep$data,fit=fits,perm=permResult,call=match.call())
        return(res)
    }else{
        return("No intervals found")
    }
}


#' @name plotTimeSeries
#' @title Plot difference function for particular bacteria
#' 
#' @details Plot the difference in abundance for 
#' 
#' @param res Output of fitTimeSeries function
#' @param C Value for which difference function has to be larger than (default 0).
#' @param xlab X-label.
#' @param ylab Y-label.
#' @param main Main label.
#' @param ... Extra plotting arguments.
#' @return NULL
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

# load("~/Dropbox/Projects/metastats/package/git/metagenomeSeq/data/mouseData.rda")
# classMatrix = aggregateByTaxonomy(mouseData,lvl='class',norm=TRUE,out='MRexperiment')
# data(mouseData)
# fitTimeSeries(obj=mouseData,feature="Actinobacteria",class="status",id="mouseID",time="relativeTime",lvl='class',B=10)
