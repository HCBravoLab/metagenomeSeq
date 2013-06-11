MRcoefs<-function(obj,by=2,coef=NULL,number=10,taxa=obj$taxa,uniqueNames=FALSE,adjust.method="fdr",group=0,eff=0,output=NULL){
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
    
    mat = cbind(tb[,coef],p)
    mat = cbind(mat,padj)
    rownames(mat) = tx;
    mat = mat[srt,]
    
    nm = c(colnames(tb)[coef],"pValue","adjPvalue")
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