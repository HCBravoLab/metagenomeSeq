MRfulltable<-function(obj,by=2,coef=NULL,number=10,taxa=obj$taxa,uniqueNames=FALSE,adjust.method="fdr",group=0,eff=0,output=NULL){
    
    tb = obj$fit$coefficients
    tx = as.character(taxa);
    
    if(uniqueNames=TRUE){
        for (nm in unique(tx)) {
            ii=which(tx==nm)
            tx[ii]=paste(tx[ii],seq_along(ii),sep=":")
        }
    }

    if(is.null(coef)){coef = 1:ncol(tb);}

    p=obj$eb$p[,by];
    padj = p.adjust(p,method=adjust.method);
    
    groups = factor(obj$fit$design[,by])
    cnts = obj$counts;
    yy = cnts>0;
    
    pa = matrix(unlist(MRfisher(obj$counts,groups,mat=TRUE)),ncol=4)
    
    np0 = rowSums(yy[,groups==unique(groups)[1]]);
    np1 = rowSums(yy[,groups==unique(groups)[2]]);

    nc0 = rowSums(cnts[,groups==unique(groups)[1]]);
    nc1 = rowSums(cnts[,groups==unique(groups)[2]]);

    if(group==0){
        srt = order(abs(tb[,by]),decreasing=TRUE)
    } else if(group==1){
        srt = order((tb[,by]),decreasing=TRUE)
    } else if(group==2){
        srt = order((tb[,by]),decreasing=FALSE)
    } else if(group==3){
        srt = order(p,decreasing=FALSE)[1:number]
    }
    valid = which(rowSums(1-obj$z)>=quantile(rowSums(1-obj$z),p=eff,na.rm=TRUE))
    srt = srt[which(srt%in%valid)][1:number]

    mat = cbind(np0,np1)
    mat = cbind(mat,nc0)
    mat = cbind(mat,nc1)
    mat = cbind(mat,pa)
    mat = cbind(mat,tb[,coef])
    mat = cbind(mat,p)
    mat = cbind(mat,padj)
    rownames(mat) = tx;
    mat = mat[srt,]

    nm = c(paste("+samples in group",unique(groups)[1]),paste("+samples in group",unique(groups)[2]),
    paste("counts in group",unique(groups)[1]),paste("counts in group",unique(groups)[2]),c("fisher.p","oddsRatio","lower","upper"),
    colnames(tb)[coef],"pValue","adjPvalue")
    colnames(mat) = nm

    if(!is.null(output)){
        nm = c("Taxa",nm)
        mat2 = cbind(rownames(mat),mat)
        mat2 = rbind(nm,mat2)
        write(t(mat2),ncolumns=ncol(mat2),file=output,sep="\t")
    } else{
        return(as.data.frame(mat))
    }
}