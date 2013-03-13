plotMRheatmap <- function(obj,n,trials,log=TRUE,norm=TRUE,...) {
  
  if(log==TRUE){
    if(norm==TRUE){
        mat = log2(cumNormMat(obj)+1);
    }else{
        mat = log2(MRcounts(obj)+1)
    }
  } else{
    if(norm==TRUE){
        mat = cumNormMat(obj)    
    }else{
        mat = MRcounts(obj)
    }
  }
  otusToKeep <- which(rowSums(mat)>0);
  otuVars=rowSds(mat[otusToKeep,]);
	otuIndices=otusToKeep[order(otuVars,decreasing=TRUE)[1:n]];
	mat2=mat[otuIndices,];
  heatmapCols=colorRampPalette(brewer.pal(9,"RdBu"))(50);
  heatmapColColors=brewer.pal(12,"Set3")[as.integer(factor(trials))];
	heatmap.2(mat2,col=heatmapCols,ColSideColors=heatmapColColors,...);
  invisible()
}