#' Table of top microbial marker gene from linear model fit including sequence
#' information
#' 
#' Extract a table of the top-ranked features from a linear model fit. This
#' function will be updated soon to provide better flexibility similar to
#' limma's topTable. This function differs from \code{link{MRcoefs}} in that it
#' provides other information about the presence or absence of features to help
#' ensure significant features called are moderately present.
#' 
#' 
#' @param obj A list containing the linear model fit produced by lmFit through
#' fitZig.
#' @param by Column number or column name specifying which coefficient or
#' contrast of the linear model is of interest.
#' @param coef Column number(s) or column name(s) specifying which coefficient
#' or contrast of the linear model to display.
#' @param number The number of bacterial features to pick out.
#' @param taxa Taxa list.
#' @param uniqueNames Number the various taxa.
#' @param adjustMethod Method to adjust p-values by. Default is "FDR". Options
#' include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#' "none". See \code{\link{p.adjust}} for more details.
#' @param group One of five choices: 0,1,2,3,4. 0: the sort is ordered by a
#' decreasing absolute value coefficient fit. 1: the sort is ordered by the raw
#' coefficient fit in decreasing order. 2: the sort is ordered by the raw
#' coefficient fit in increasing order. 3: the sort is ordered by the p-value
#' of the coefficient fit in increasing order. 4: no sorting.
#' @param eff Filter features to have at least a "eff" quantile or number of effective samples.
#' @param numberEff Boolean, whether eff should represent quantile (default/FALSE) or number.
#' @param ncounts Filter features to those with at least 'counts' counts.
#' @param file Name of output file, including location, to save the table.
#' @return Table of the top-ranked features determined by the linear fit's
#' coefficient.
#' @seealso \code{\link{fitZig}} \code{\link{fitFeatureModel}} \code{\link{MRcoefs}} \code{\link{MRtable}}
#' \code{\link{fitPA}}
#' @examples
#' 
#' data(lungData)
#' k = grep("Extraction.Control",pData(lungData)$SampleType)
#' lungTrim = lungData[,-k]
#' lungTrim=filterData(lungTrim,present=30)
#' lungTrim=cumNorm(lungTrim,p=0.5)
#' smokingStatus = pData(lungTrim)$SmokingStatus
#' mod = model.matrix(~smokingStatus)
#' fit = fitZig(obj = lungTrim,mod=mod)
#' # head(MRfulltable(fit))
#' ####
#' fit = fitFeatureModel(obj = lungTrim,mod=mod)
#' # head(MRfulltable(fit))
#'
MRfulltable<-function(obj,by=2,coef=NULL,number=10,taxa=obj$taxa,
    uniqueNames=FALSE,adjustMethod="fdr",group=0,eff=0,numberEff=FALSE,ncounts=0,file=NULL){
    
    if(length(grep("fitFeatureModel",obj$call))){
        groups = factor(obj$design[,by])
        by = "logFC"; coef = 1:2;
        tb = data.frame(logFC=obj$fitZeroLogNormal$logFC,se=obj$fitZeroLogNormal$se)
        p  = obj$pvalues
    } else {
        tb = obj$fit$coefficients
        if(is.null(coef)){
            coef = 1:ncol(tb)
        }
        p=obj$eb$p.value[,by]
        groups = factor(obj$fit$design[,by])
        if(eff>0){
            effectiveSamples = calculateEffectiveSamples(obj)
            if(numberEff == FALSE){
                valid = which(effectiveSamples>=quantile(effectiveSamples,p=eff,na.rm=TRUE))
            } else {
                valid = which(effectiveSamples>=eff)
            }
        }
    }

    tx = as.character(taxa)
    if(uniqueNames==TRUE){
        for (nm in unique(tx)) {
            ii=which(tx==nm)
            tx[ii]=paste(tx[ii],seq_along(ii),sep=":")
        }
    }
    padj = p.adjust(p,method=adjustMethod)
    cnts = obj$counts
    yy = cnts>0
    
    pa = matrix(unlist(fitPA(obj$counts,groups)),ncol=5)
    
    np0 = rowSums(yy[,groups==0])
    np1 = rowSums(yy[,groups==1])

    nc0 = rowSums(cnts[,groups==0])
    nc1 = rowSums(cnts[,groups==1])

    if(group==0){
        srt = order(abs(tb[,by]),decreasing=TRUE)
    } else if(group==1){
        srt = order((tb[,by]),decreasing=TRUE)
    } else if(group==2){
        srt = order((tb[,by]),decreasing=FALSE)
    } else if(group==3){
        srt = order(p,decreasing=FALSE)
    } else {
        srt = 1:length(padj)
    }

    valid = 1:length(padj)
    if(ncounts>0){
        np=rowSums(cbind(np0,np1))
        valid = intersect(valid,which(np>=ncounts))
    }
    srt = srt[which(srt%in%valid)][1:min(number,nrow(tb))]

    mat = cbind(np0,np1)
    mat = cbind(mat,nc0)
    mat = cbind(mat,nc1)
    mat = cbind(mat,pa)
    mat = cbind(mat,tb[,coef])
    mat = cbind(mat,p)
    mat = cbind(mat,padj)
    rownames(mat) = tx
    mat = mat[srt,]

    nm = c("+samples in group 0","+samples in group 1","counts in group 0",
        "counts in group 1",c("oddsRatio","lower","upper","fisherP","fisherAdjP"),
    colnames(tb)[coef],"pvalues","adjPvalues")
    colnames(mat) = nm

    if(!is.null(file)){
        nm = c("Taxa",nm)
        mat2 = cbind(rownames(mat),mat)
        mat2 = rbind(nm,mat2)
        write(t(mat2),ncolumns=ncol(mat2),file=file,sep="\t")
    }
    return(as.data.frame(mat))
}
