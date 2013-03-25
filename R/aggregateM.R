#' Aggregates the counts to a particular classification.
#'
#' This function takes an eSet object of data at a particular level with feature information allowing
#' for aggregation of counts to a particular level. This method assumes taxa begin at the highest level and continue to the current level.
#'
#' @param obj An eSet object of count data.
#' @param lvl The level to go up (numeric, 1,2,3).
#' @param split The way character strings in taxa in the obj are split.
#' @return Updated eSet object with counts aggregated to the various taxanomic levels.
#'
#' @name aggregateM
aggregateM <-
function(obj,taxa,lvl,split=";"){

	tmp<-strsplit(taxa,split=split);
	nrows = length(taxa);
    cnts  = MRcounts(obj)
	
	maxName=max(sapply(tmp,length))
	
	featureMat<-array("NA",dim=c(nrows,maxName));
	
	for(i in 1:nrows){
		for(j in 1:length(tmp[[i]])){
			featureMat[i,j]=tmp[[i]][j]
			}
	}
	
	tt = sapply(1:nrows,function(i){
		   t = featureMat[i,1]
		   for(j in 2:lvl){
				t = paste(t,featureMat[i,j],sep=";")
		   }
		   t})
	
	newTaxa<-unique(tt);
    newCountMat<-array(0,dim=c(length(newTaxa),ncol(cnts)));

    colnames(newCountMat)<-sampleNames(obj)
	rownames(newCountMat)<-newTaxa;	
	
	for(i in 1:length(newTaxa)){
		if(length(which(tt==newTaxa[i]))<2){
			newCountMat[rownames(newCountMat)==newTaxa[i],]=as.matrix(cnts[which(tt==newTaxa[i]),],nr=1)
		}else{
			newCountMat[rownames(newCountMat)==newTaxa[i],]=colSums(cnts[which(tt==newTaxa[i]),])
		}
	}
	
	taxa<-newTaxa;
	assayData(obj)<-newCountMat;
	validObject(obj)

	return(obj)
}
