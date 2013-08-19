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
#' @return NA
#' @seealso \code{\link{cumNorm}} \code{\link{fitZig}}
#' @examples
#' 
#' data(lungData)
#' k = grep("Extraction.Control",pData(lungData)$SampleType)
#' lungTrim = lungData[,-k]
#' lungTrim = lungTrim[-which(rowSums(MRcounts(lungTrim)>0)<20),]
#' res = MRfisher(lungTrim,pData(lungTrim)$SmokingStatus);
#' head(res)
#' 
MRfisher<-function(obj,cl,thres=0){
    if(class(obj)=="MRexperiment"){
        x = MRcounts(obj)>thres;
    } else if(class(obj) == "matrix") {
        x = obj>thres;
    } else {
        stop("Object needs to be either a MRexperiment object or matrix")
    }
    nrows= nrow(x);
	if(is.null(rownames(x))){rownames(x)=1:nrows}
    
    res = sapply(1:nrows,function(i){
        tbl = table(1-x[i,],cl)
        if(sum(dim(tbl))!=4){
            tbl = array(0,dim=c(2,2));
            tbl[1,1] = sum(x[i,cl==unique(cl)[1]])
            tbl[1,2] = sum(x[i,cl==unique(cl)[2]])
            tbl[2,1] = sum(cl==unique(cl)[1])-tbl[1,1]
            tbl[2,2] = sum(cl==unique(cl)[2])-tbl[1,2]
        }
        ft <- fisher.test(tbl, workspace = 8e6, alternative = "two.sided", conf.int = T)
        list(p=ft$p.value,o=ft$estimate,cl=ft$conf.int[1],cu=ft$conf.int[2])
    })
    
    dat = data.frame(as.matrix(t(res)))
    rownames(dat) = rownames(x)
    colnames(dat) = c("pvalues","oddsRatio","lower","upper")
    return(dat)
}
