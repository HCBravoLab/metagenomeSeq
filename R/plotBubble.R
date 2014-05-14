#' Basic plot of binned vectors.
#' 
#' This function plots takes two vectors, calculates the contingency table and 
#' plots circles sized by the contingency table value. Optional significance vectors
#' of the values significant will shade the circles by proportion of significance.
#' 
#' 
#' @param yvector A vector of values represented along y-axis.
#' @param xvector A vector of values represented along x-axis.
#' @param sigvector A vector of the names of significant features (names should match x/yvector).
#' @param nbreaks Number of bins to break yvector and xvector into.
#' @param ret Boolean to return the observed data that would have been plotted.
#' @param scale Scaling of circle bins.
#' @param ... Additional plot arguments.
#' @return NA
#' @seealso \code{\link{plotMRheatmap}}
#' @examples
#' 
#' data(mouseData)
#' mouseData = mouseData[which(rowSums(mouseData)>139),]
#' sparsity = rowMeans(MRcounts(mouseData)==0)
#' lor = log(fitPA(mouseData,cl=pData(mouseData)[,3])$oddsRatio)
#' plotBubble(lor,sparsity,main="lor ~ sparsity")
#' 
plotBubble<-function(yvector,xvector,sigvector=NULL,nbreaks=10,ret=FALSE,scale=1,...){
    ybreaks = cut(yvector,breaks=quantile(yvector,p=seq(0,1,length.out=nbreaks)),include.lowest=T)
    xbreaks = cut(xvector,breaks=quantile(xvector,p=seq(0,1,length.out=nbreaks)),include.lowest=T)
    contTable = lapply(levels(xbreaks),function(i){
        k = which(xbreaks==i)
        sapply(levels(ybreaks),function(j){
            length(which(ybreaks[k]==j))
            })
        })
    names(contTable) = levels(xbreaks)
    yvec = 1:length(levels(ybreaks))
    nc = length(yvec)

    if(!is.null(sigvector)){
        # I am calculating contTable twice if sigvector==TRUE
        # This can be changed to if else statement to return two rows
        contSig = lapply(levels(xbreaks),function(i){
            k = which(xbreaks==i)
            sapply(levels(ybreaks),function(j){
                x = sum(names(yvector[k])[which(ybreaks[k]==j)]%in%sigvector)/length(which(ybreaks[k]==j))
                if(is.na(x)) x = 0
                x
            })
        })
    } else {
        contSig = lapply(levels(xbreaks),function(i){
            k = which(xbreaks==i)
            sapply(levels(ybreaks),function(j){
                1
            })
        })
    }
    
    medianSizes = median(unlist(contTable))
    plot(y=yvec,x=rep(1,nc),cex=scale*contTable[[1]]/medianSizes,
        xlim=c(-0.25,nc+.25),ylim=c(-0.25,nc+.25),bty="n",xaxt="n",yaxt="n",
        xlab="",ylab="",pch=21,...,bg=rgb(blue=1,red=0,green=0,alpha=contSig[[1]]))
    for(i in 2:length(contTable)){
        points(y=yvec,x=rep(i,nc),cex =scale*contTable[[i]]/medianSizes,pch=21,bg=rgb(blue=1,red=0,green=0,alpha=contSig[[i]]))
    }
    axis(1,at = 1:nc,labels=levels(xbreaks),las=2,cex.axis=.5)
    axis(2,at = 1:nc,labels=levels(ybreaks),las=2,cex.axis=.5)
    if(ret == TRUE){
        return(cbind(xbreaks,ybreaks))
    }
}
