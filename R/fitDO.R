#' Wrapper to calculate Discovery Odds Ratios on feature values.
#' 
#' This function returns a data frame of p-values, odds ratios, lower and upper
#' confidence limits for every row of a matrix. The discovery odds ratio is calculated
#' as using Fisher's exact test on actual counts. The test's hypothesis is whether 
#' or not the discovery of counts for a feature (of all counts) is found in greater proportion
#' in a particular group.
#' 
#' 
#' @param obj A MRexperiment object with a count matrix, or a simple count
#' matrix.
#' @param cl Group comparison
#' @param norm Whether or not to normalize the counts - if MRexperiment object.
#' @param log Whether or not to log2 transform the counts - if MRexperiment object.
#' @param parallel Use multiple cores?
#' @param cores Number of cores to use.
#' @param adjust.method Method to adjust p-values by. Default is "FDR". Options
#' include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",                                                                                                                            
#' "none". See \code{\link{p.adjust}} for more details.
#' @param ... Extra options for makeCluster
#' @return Matrix of odds ratios, p-values, lower and upper confidence intervals
#' @seealso \code{\link{cumNorm}} \code{\link{fitZig}} \code{\link{fitPA}} \code{\link{fitMeta}}
#' @examples
#' 
#' data(lungData)
#' k = grep("Extraction.Control",pData(lungData)$SampleType)
#' lungTrim = lungData[,-k]
#' lungTrim = lungTrim[-which(rowSums(MRcounts(lungTrim)>0)<20),]
#' res = fitDO(lungTrim,pData(lungTrim)$SmokingStatus);
#' head(res)
#' 
fitDO<-function(obj,cl,norm=TRUE,log=TRUE,parallel=FALSE,cores=2,adjust.method='fdr',...){
    if(class(obj)=="MRexperiment"){
        x = MRcounts(obj,norm=norm,log=log);
    } else if(class(obj) == "matrix") {
        x = obj;
    } else {
        stop("Object needs to be either a MRexperiment object or matrix")
    }
    nrows= nrow(x);
	if(is.null(rownames(x))){rownames(x)=1:nrows}

    sumClass1 = round(sum(x[,cl==unique(cl)[1]]))
    sumClass2 = round(sum(x[,cl==unique(cl)[2]]))
    if(parallel==FALSE){
        res = sapply(1:nrows,function(i){
            tbl = table(1-x[i,],cl)
            if(sum(dim(tbl))!=4){
                tbl = array(0,dim=c(2,2));
                tbl[1,1] = round(sum(x[i,cl==unique(cl)[1]]))
                tbl[1,2] = round(sum(x[i,cl==unique(cl)[2]]))
                tbl[2,1] = sumClass1-tbl[1,1]
                tbl[2,2] = sumClass2-tbl[1,2]
            }
            ft <- fisher.test(tbl, workspace = 8e6, alternative = "two.sided", conf.int = TRUE)
            cbind(o=ft$estimate,cl=ft$conf.int[1],cu=ft$conf.int[2],p=ft$p.value)
        })
        res = data.frame(as.matrix(t(res)))
        adjp = p.adjust(res[,4],method=adjust.method)
        res = data.frame(cbind(res,adjp))
    } else {
        library(parallel)
        cores <- makeCluster(getOption("cl.cores", cores))
        res = parRapply(cl=cores,x,function(i){
                tbl = table(1-i,cl)
                if(sum(dim(tbl))!=4){
                    tbl = array(0,dim=c(2,2));
                    tbl[1,1] = round(sum(i[cl==unique(cl)[1]]))
                    tbl[1,2] = round(sum(i[cl==unique(cl)[2]]))
                    tbl[2,1] = sumClass1-tbl[1,1]
                    tbl[2,2] = sumClass2-tbl[1,2]
                }
                ft <- fisher.test(tbl,workspace=8e6,alternative="two.sided",conf.int=TRUE)
                cbind(p=ft$p.value,o=ft$estimate,cl=ft$conf.int[1],cu=ft$conf.int[2])
            })
        stopCluster(cores)
        nres = nrows*4
        seqs = seq(1,nres,by=4)
        p = res[seqs]
        adjp = p.adjust(p,method=adjust.method)
        o = res[seqs+1]
        cl = res[seqs+2]
        cu = res[seqs+3]
        res = data.frame(cbind(o,cl,cu,p,adjp))
    }
    colnames(res) = c("oddsRatio","lower","upper","pvalues","adjPvalues")
    rownames(res) = rownames(x)
    return(res)
}
