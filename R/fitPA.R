#' Wrapper to run fisher's test on presence/absence of a feature.
#' 
#' This function returns a data frame of p-values, odds ratios, lower and upper
#' confidence limits for every row of a matrix.
#' 
#' 
#' @param obj A MRexperiment object with a count matrix, or a simple count
#' matrix.
#' @param cl Group comparison
#' @param thres Threshold for defining presence/absence.
#' @param parallel Use multiple cores?
#' @param cores Number of cores to use.
#' @param adjust.method Method to adjust p-values by. Default is "FDR". Options
#' include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#' "none". See \code{\link{p.adjust}} for more details.
#' @param ... Extra options for makeCluster
#' @return Matrix of odds ratios, p-values, lower and upper confidence intervals
#' @seealso \code{\link{cumNorm}} \code{\link{fitZig}} \code{\link{fitDO}} \code{\link{fitMeta}}
#' @examples
#' 
#' data(lungData)
#' k = grep("Extraction.Control",pData(lungData)$SampleType)
#' lungTrim = lungData[,-k]
#' lungTrim = lungTrim[-which(rowSums(MRcounts(lungTrim)>0)<20),]
#' res = fitPA(lungTrim,pData(lungTrim)$SmokingStatus);
#' head(res)
#'
fitPA<-function(obj,cl,thres=0,parallel=FALSE,cores=2,adjust.method='fdr',...){
    if(class(obj)=="MRexperiment"){
        x = MRcounts(obj)>thres;
    } else if(class(obj) == "matrix") {
        x = obj>thres;
    } else {
        stop("Object needs to be either a MRexperiment object or matrix")
    }
    nrows= nrow(x);
    if(is.null(rownames(x))){rownames(x)=1:nrows}

    nClass1 = sum(cl==unique(cl)[1])
    nClass2 = sum(cl==unique(cl)[2])

    if(parallel==FALSE){
        res = sapply(1:nrows,function(i){
            tbl = table(1-x[i,],cl)
            if(sum(dim(tbl))!=4){
                tbl = array(0,dim=c(2,2));
                tbl[1,1] = sum(x[i,cl==unique(cl)[1]])
                tbl[1,2] = sum(x[i,cl==unique(cl)[2]])
                tbl[2,1] = nClass1-tbl[1,1]
                tbl[2,2] = nClass2-tbl[1,2]
            }
            ft <- fisher.test(tbl,workspace=8e6,alternative="two.sided",conf.int=TRUE)
            cbind(o=ft$estimate,cl=ft$conf.int[1],cu=ft$conf.int[2],p=ft$p.value)
        })
        res = data.frame(as.matrix(t(res)))
        adjp = p.adjust(res[,"p"],method=adjust.method)
        res = data.frame(cbind(res,adjp))
    } else {
        library(parallel)
        cores <- makeCluster(getOption("cl.cores", cores))
        res = parRapply(cl=cores,x,function(i){
                tbl = table(1-i,cl)
                if(sum(dim(tbl))!=4){
                    tbl = array(0,dim=c(2,2));
                    tbl[1,1] = sum(i[cl==unique(cl)[1]])
                    tbl[1,2] = sum(i[cl==unique(cl)[2]])
                    tbl[2,1] = nClass1-tbl[1,1]
                    tbl[2,2] = nClass2-tbl[1,2]
                }
                ft <- fisher.test(tbl,workspace=8e6,alternative="two.sided",conf.int=TRUE)
                cbind(o=ft$estimate,cl=ft$conf.int[1],cu=ft$conf.int[2],p=ft$p.value)
            })
        stopCluster(cores)
        nres = nrows*4
        seqs = seq(1,nres,by=4)
        p = res[seqs+3]
        adjp = p.adjust(p,method=adjust.method)
        o = res[seqs]
        cl = res[seqs+1]
        cu = res[seqs+2]
        res = data.frame(cbind(o,cl,cu,p,adjp))
    }
    colnames(res) = c("oddsRatio","lower","upper","pvalues","adjPvalues")
    rownames(res) = rownames(x)
    return(res)
}
