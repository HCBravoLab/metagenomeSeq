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
#' @return NA
#' @seealso \code{\link{cumNorm}} \code{\link{fitZig}}
#' @examples
#' 
#' data(lungData)
#' k = grep("Extraction.Control",pData(lungData)$SampleType)
#' lungTrim = lungData[,-k]
#' lungTrim = lungTrim[-which(rowSums(MRcounts(lungTrim)>0)<20),]
#' res = fitDO(lungTrim,pData(lungTrim)$SmokingStatus);
#' head(res)
#' 
fitDO<-function(obj,cl,norm=TRUE,log=TRUE){
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

    res = sapply(1:nrows,function(i){
        tbl = table(1-x[i,],cl)
        if(sum(dim(tbl))!=4){
            tbl = array(0,dim=c(2,2));
            tbl[1,1] = round(sum(x[i,cl==unique(cl)[1]]))
            tbl[1,2] = round(sum(x[i,cl==unique(cl)[2]]))
            tbl[2,1] = sumClass1-tbl[1,1]
            tbl[2,2] = sumClass2-tbl[1,2]
        }
        ft <- fisher.test(tbl, workspace = 8e6, alternative = "two.sided", conf.int = T)
        list(p=ft$p.value,o=ft$estimate,cl=ft$conf.int[1],cu=ft$conf.int[2])
    })
    
    dat = data.frame(as.matrix(t(res)))
    rownames(dat) = rownames(x)
    colnames(dat) = c("pvalues","oddsRatio","lower","upper")
    return(dat)
}
