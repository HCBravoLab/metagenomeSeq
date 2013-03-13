pkgname <- "metagenomeSeq"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('metagenomeSeq')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("MRcoefs")
### * MRcoefs

flush(stderr()); flush(stdout())

### Name: MRcoefs
### Title: Table of top-ranked microbial marker gene from linear model fit
### Aliases: MRcoefs

### ** Examples

data(lungData)
k = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-k]
k = which(rowSums(MRcounts(lungTrim)>0)<10)
lungTrim = lungTrim[-k,]
cumNorm(lungTrim)
smokingStatus = pData(lungTrim)$SmokingStatus
mod = model.matrix(~smokingStatus)
settings = zigControl(maxit=1,verbose=FALSE)
fit = fitZig(obj = lungTrim,mod=mod,control=settings)
head(MRcoefs(fit))



cleanEx()
nameEx("MRcounts")
### * MRcounts

flush(stderr()); flush(stdout())

### Name: MRcounts
### Title: Accessor for the counts slot of a MRexperiment object
### Aliases: MRcounts,MRexperiment-method MRcounts

### ** Examples

   data(lungData)
   head(MRcounts(lungData))



cleanEx()
nameEx("MRexperiment-class")
### * MRexperiment-class

flush(stderr()); flush(stdout())

### Name: MRexperiment
### Title: Class "MRexperiment" - a modified eSet object for the data from
###   high-throughput sequencing experiments
### Aliases: MRexperiment-class [,MRexperiment-method

### ** Examples

# See vignette



cleanEx()
nameEx("MRfisher")
### * MRfisher

flush(stderr()); flush(stdout())

### Name: MRfisher
### Title: Wrapper to run fisher's test on presence/absence of a feature.
### Aliases: MRfisher

### ** Examples

data(lungData)
k = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-k]
lungTrim = lungTrim[-which(rowSums(MRcounts(lungTrim)>0)<20),]
res = MRfisher(lungTrim,pData(lungTrim)$SmokingStatus);
head(res)



cleanEx()
nameEx("MRfulltable")
### * MRfulltable

flush(stderr()); flush(stdout())

### Name: MRfulltable
### Title: Table of top microbial marker gene from linear model fit
###   including sequence information
### Aliases: MRfulltable

### ** Examples

data(lungData)
k = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-k]
k = which(rowSums(MRcounts(lungTrim)>0)<10)
lungTrim = lungTrim[-k,]
cumNorm(lungTrim)
smokingStatus = pData(lungTrim)$SmokingStatus
mod = model.matrix(~smokingStatus)
settings = zigControl(maxit=1,verbose=FALSE)
fit = fitZig(obj = lungTrim,mod=mod,control=settings)
head(MRfulltable(fit))



cleanEx()
nameEx("MRtable")
### * MRtable

flush(stderr()); flush(stdout())

### Name: MRtable
### Title: Table of top microbial marker gene from linear model fit
###   including sequence information
### Aliases: MRtable

### ** Examples

data(lungData)
k = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-k]
k = which(rowSums(MRcounts(lungTrim)>0)<10)
lungTrim = lungTrim[-k,]
cumNorm(lungTrim)
smokingStatus = pData(lungTrim)$SmokingStatus
mod = model.matrix(~smokingStatus)
settings = zigControl(maxit=1,verbose=FALSE)
fit = fitZig(obj = lungTrim,mod=mod,control=settings)
head(MRtable(fit))



cleanEx()
nameEx("cumNorm")
### * cumNorm

flush(stderr()); flush(stdout())

### Name: cumNorm
### Title: Cumulative sum scaling factors.
### Aliases: cumNorm

### ** Examples

data(mouseData)
cumNorm(mouseData)
head(normFactors(mouseData))



cleanEx()
nameEx("cumNormMat")
### * cumNormMat

flush(stderr()); flush(stdout())

### Name: cumNormMat
### Title: Cumulative sum scaling factors.
### Aliases: cumNormMat

### ** Examples

data(mouseData)
head(cumNormMat(mouseData))



cleanEx()
nameEx("cumNormStat")
### * cumNormStat

flush(stderr()); flush(stdout())

### Name: cumNormStat
### Title: Cumulative sum scaling percentile selection
### Aliases: cumNormStat

### ** Examples

data(mouseData)
p = round(cumNormStat(mouseData,pFlag=FALSE),digits=2)
s95=cumNorm(mouseData)



cleanEx()
nameEx("expSummary")
### * expSummary

flush(stderr()); flush(stdout())

### Name: expSummary
### Title: Access MRexperiment object experiment data
### Aliases: expSummary,MRexperiment-method expSummary

### ** Examples

data(mouseData)
expSummary(mouseData)



cleanEx()
nameEx("exportMat")
### * exportMat

flush(stderr()); flush(stdout())

### Name: exportMat
### Title: export the normalized eSet dataset as a matrix.
### Aliases: exportMatrix exportMat

### ** Examples

# see vignette



cleanEx()
nameEx("exportStats")
### * exportStats

flush(stderr()); flush(stdout())

### Name: exportStats
### Title: Various statistics of the count data.
### Aliases: exportStats

### ** Examples

# see vignette



cleanEx()
nameEx("fitZig")
### * fitZig

flush(stderr()); flush(stdout())

### Name: fitZig
### Title: Computes the weighted fold-change estimates and t-statistics.
### Aliases: fitZig

### ** Examples

data(lungData)
k = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-k]
cumNorm(lungTrim)
k = which(rowSums(MRcounts(lungTrim)>0)<10)
lungTrim = lungTrim[-k,]
smokingStatus = pData(lungTrim)$SmokingStatus
mod = model.matrix(~smokingStatus)
settings = zigControl(maxit=1,verbose=FALSE)
fit = fitZig(obj = lungTrim,mod=mod,control=settings)



cleanEx()
nameEx("libSize")
### * libSize

flush(stderr()); flush(stdout())

### Name: libSize
### Title: Access sample depth of coverage from MRexperiment object
### Aliases: libSize,MRexperiment-method libSize

### ** Examples

   data(lungData)
   head(libSize(lungData))



cleanEx()
nameEx("load_meta")
### * load_meta

flush(stderr()); flush(stdout())

### Name: load_meta
### Title: Load a count dataset associated with a study.
### Aliases: load_meta metagenomicLoader

### ** Examples

dataDirectory <- system.file("extdata", package="metagenomeSeq")
lung = load_meta(file.path(dataDirectory,"CHK_NAME.otus.count.csv"))



cleanEx()
nameEx("load_metaQ")
### * load_metaQ

flush(stderr()); flush(stdout())

### Name: load_metaQ
### Title: Load a count dataset associated with a study set up in a Qiime
###   format.
### Aliases: load_metaQ qiimeLoader

### ** Examples

# see vignette



cleanEx()
nameEx("load_phenoData")
### * load_phenoData

flush(stderr()); flush(stdout())

### Name: load_phenoData
### Title: Load a clinical/phenotypic dataset associated with a study.
### Aliases: load_phenoData phenoData

### ** Examples

# see vignette



cleanEx()
nameEx("newMRexperiment")
### * newMRexperiment

flush(stderr()); flush(stdout())

### Name: newMRexperiment
### Title: Create a MRexperiment object
### Aliases: newMRexperiment

### ** Examples

cnts = matrix(abs(rnorm(1000)),nc=10)
obj <- newMRexperiment(cnts)



cleanEx()
nameEx("normFactors")
### * normFactors

flush(stderr()); flush(stdout())

### Name: normFactors
### Title: Access the normalization factors in a MRexperiment object
### Aliases: normFactors,MRexperiment-method normFactors

### ** Examples

   data(lungData)
   cumNorm(lungData)
   head(normFactors(lungData))



cleanEx()
nameEx("plotCorr")
### * plotCorr

flush(stderr()); flush(stdout())

### Name: plotCorr
### Title: Basic correlation plot function for normalized or unnormalized
###   counts.
### Aliases: plotCorr

### ** Examples

data(mouseData)
trials = pData(mouseData)$diet
plotCorr(obj=mouseData,n=200,cexRow = 0.4,cexCol = 0.4,trace="none",dendrogram="none")



cleanEx()
nameEx("plotGenus")
### * plotGenus

flush(stderr()); flush(stdout())

### Name: plotGenus
### Title: Basic plot function of the raw or normalized data.
### Aliases: genusPlot plotGenus

### ** Examples

data(mouseData)
classIndex=list(controls=which(pData(mouseData)$diet=="BK"))
classIndex$cases=which(pData(mouseData)$diet=="Western")
otuIndex = grep("Strep",fData(mouseData)$fdata)
otuIndex=otuIndex[order(rowSums(MRcounts(mouseData)[otuIndex,]),decreasing=TRUE)]
plotGenus(mouseData,otuIndex,classIndex,xlab="OTU log-normalized counts",no=1:2,xaxt="n",norm=FALSE,ylab="Strep")
lablist<-rep(c("Controls","Cases"),times=2)
axis(1, at=seq(1,4,by=1), labels = lablist)




cleanEx()
nameEx("plotMRheatmap")
### * plotMRheatmap

flush(stderr()); flush(stdout())

### Name: plotMRheatmap
### Title: Basic heatmap plot function for normalized counts.
### Aliases: plotMRheatmap

### ** Examples

data(mouseData)
trials = pData(mouseData)$diet
plotMRheatmap(obj=mouseData,n=200,trials=trials,cexRow = 0.4,cexCol = 0.4,trace="none")



cleanEx()
nameEx("plotOTU")
### * plotOTU

flush(stderr()); flush(stdout())

### Name: plotOTU
### Title: Basic plot function of the raw or normalized data.
### Aliases: plotOTU

### ** Examples

data(mouseData)
classIndex=list(controls=which(pData(mouseData)$diet=="BK"))
classIndex$cases=which(pData(mouseData)$diet=="Western")
# you can specify whether or not to normalize, and to what level
plotOTU(mouseData,otu=9083,classIndex,xlab="OTU log-normalized counts",norm=FALSE,xaxt="n",main="9083 feature abundances")
lablist<- c("Controls","Cases")
axis(1, at=seq(1,2,by=1), labels = lablist)



cleanEx()
nameEx("posterior.probs")
### * posterior.probs

flush(stderr()); flush(stdout())

### Name: posterior.probs
### Title: Access the posterior probabilities that results from analysis
### Aliases: posterior.probs,MRexperiment-method posterior.probs

### ** Examples

# see vignette



cleanEx()
nameEx("zigControl")
### * zigControl

flush(stderr()); flush(stdout())

### Name: zigControl
### Title: Settings for the fitZig function
### Aliases: settings zigControl

### ** Examples

control =  zigControl(tol=1e-10,maxit=10,verbose=FALSE)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
