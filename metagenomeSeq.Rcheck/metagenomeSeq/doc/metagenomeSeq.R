### R code from vignette source 'metagenomeSeq.Rnw'

###################################################
### code chunk number 1: config
###################################################
options(width = 50)
options(continue=" ")
options(prompt="R> ")
set.seed(42)


###################################################
### code chunk number 2: loadMatrix
###################################################
library(metagenomeSeq)
dataDirectory <- system.file("extdata", package="metagenomeSeq")
lung = load_meta(file.path(dataDirectory,"CHK_NAME.otus.count.csv")) 
dim(lung$counts)


###################################################
### code chunk number 3: loadTaxa
###################################################
taxa = read.csv(file.path(dataDirectory,"CHK_otus.taxonomy.csv"),
                sep="\t",header=T,stringsAsFactors=F)[,2]
otu  = read.csv(file.path(dataDirectory,"CHK_otus.taxonomy.csv"),
                sep="\t",header=T,stringsAsFactors=F)[,1]


###################################################
### code chunk number 4: loadClin
###################################################
clin = load_phenoData(file.path(dataDirectory,"CHK_clinical.csv"),tran=TRUE)
ord = match(colnames(lung$counts),rownames(clin)) 
clin = clin[ord,]
head(clin[1:2,])


###################################################
### code chunk number 5: createMRexperiment1
###################################################
phenotypeData = as(clin,"AnnotatedDataFrame")
phenotypeData


###################################################
### code chunk number 6: createMRexperiment2
###################################################
OTUdata = as(lung$taxa,"AnnotatedDataFrame")
varLabels(OTUdata) = "taxa"
OTUdata


###################################################
### code chunk number 7: createMRexperiment2
###################################################
counts = lung$counts
obj = newMRexperiment(counts,phenoData=phenotypeData,featureData=OTUdata)
#experimentData(obj) = annotate::pmid2MIAME("21680950")
obj


###################################################
### code chunk number 8: calculateP
###################################################
data(lungData)
p=cumNormStat(lungData)


###################################################
### code chunk number 9: normalizeData
###################################################
cumNorm(lungData,p=p)


###################################################
### code chunk number 10: saveData
###################################################
mat = MRcounts(lungData,norm=TRUE)[1:5,1:5]
exportMat(mat,output=file.path(dataDirectory,"temp.tsv"))


###################################################
### code chunk number 11: saveData2
###################################################
exportStats(lungData[,1:5],output=file.path(dataDirectory,"temp.tsv"),p=p)
head(read.csv(file=file.path(dataDirectory,"temp.tsv"),sep="\t"))


###################################################
### code chunk number 12: removeData
###################################################
system(paste("rm",file.path(dataDirectory,"temp.tsv")))


###################################################
### code chunk number 13: preprocess
###################################################
data(lungData)
k = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-k]
k = which(rowSums(MRcounts(lungTrim)>0)<10)
lungTrim = lungTrim[-k,]
cumNorm(lungTrim)


###################################################
### code chunk number 14: Model
###################################################
smokingStatus = pData(lungTrim)$SmokingStatus
bodySite = pData(lungTrim)$SampleType
mod = model.matrix(~smokingStatus+bodySite)
settings = zigControl(maxit=10,verbose=TRUE)
fit = fitZig(obj = lungTrim,mod=mod,control=settings)


###################################################
### code chunk number 15: fittedResult
###################################################
taxa = 
  sapply(strsplit(as.character(fData(lungTrim)$taxa),split=";"),
         function(i){i[length(i)]})
head(MRcoefs(fit,taxa=taxa,coef=2))


###################################################
### code chunk number 16: heatmapData
###################################################
data(mouseData)
trials = pData(mouseData)$diet
plotMRheatmap(obj=mouseData,n=200,trials=trials,
              cexRow = 0.4,cexCol = 0.4,trace="none")
plotCorr(obj=mouseData,n=200,
              cexRow = 0.25,cexCol = 0.25,trace="none",dendrogram="none")


###################################################
### code chunk number 17: plotOTUData
###################################################
head(MRtable(fit,coef=2,taxa=1:length(fData(lungTrim)$taxa)))
patients=sapply(strsplit(rownames(pData(lungTrim)),split="_"),function(i){i[3]})
pData(lungTrim)$patients=patients
classIndex=list(controls=which(pData(lungTrim)$SmokingStatus=="Smoker"))
classIndex$cases=which(pData(lungTrim)$SmokingStatus=="NonSmoker")
otu = 779
x = fData(lungTrim)$taxa[otu]
plotOTU(lungTrim,otu=otu,classIndex,xaxt="n",
        ylab="Normalized log(cpt)",main="Neisseria meningitidis")
lablist<- c("Smoker","NonSmoker")
axis(1, at=seq(1,2,by=1), labels = lablist)


###################################################
### code chunk number 18: plotGenusData
###################################################
otulist = grep(x,fData(lungTrim)$taxa)
plotGenus(lungTrim,otulist,classIndex,xaxt="n",
          ylab="Normalized log(cpt)",main="Neisseria meningitidis")
lablist<- c("S","NS")
axis(1, at=seq(1,6,by=1), labels = rep(lablist,times=3))


###################################################
### code chunk number 19: sessi
###################################################
sessionInfo()


