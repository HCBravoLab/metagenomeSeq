#' MRexperiment to biom objects
#' 
#' Wrapper to convert MRexperiment objects to biom objects.
#' 
#' @param obj The MRexperiment object.
#' @param id Optional id for the biom matrix.
#' @param norm normalize count table
#' @param log log2 transform count table
#' @param sl scaling factor for normalized counts.
#' @param qiimeVersion Format fData according to QIIME specifications (assumes only taxonomy in fData).
#' @return A biom object.
#' @seealso \code{\link{loadMeta}} \code{\link{loadPhenoData}} \code{\link{newMRexperiment}} \code{\link{loadBiom}} \code{\link{biom2MRexperiment}}
MRexperiment2biom <- function(obj,id=NULL,norm=FALSE,log=FALSE,sl=1000,qiimeVersion=TRUE){
    requireNamespace("biomformat")
    id = id
    format = "Biological Observation Matrix 1.0.0-dev"
    format_url = "http://biom-format.org/documentation/format_versions/biom-1.0.html"
    type = "OTU table"
    generated_by = sprintf("metagenomeSeq %s",packageVersion("metagenomeSeq"))
    date = as.character(Sys.time())
    matrix_type = "dense"
    matrix_element_type = "int"
    if( (norm==TRUE) | (log == TRUE) ) {
        matrix_element_type = "float"
    }
    
    data  = MRcounts(obj,norm=norm,log=log,sl=sl)
    shape = dim(data)

    rows   = metadata(fData(obj),qiimeVersion=qiimeVersion)
    columns= metadata(pData(obj))

    data = as.list(as.data.frame(t(data)))
    names(data) <- NULL

    biomlist = list(id=id,format=format,format_url=format_url,type=type,generated_by=generated_by,
                    date=date,matrix_type=matrix_type,matrix_element_type=matrix_element_type,shape=shape,
                    rows=rows,columns=columns,data=data)
    biomformat::biom(biomlist)
}

metadata <- function(df,qiimeVersion=FALSE){
    if(ncol(df)>0){
        for(i in 1:ncol(df)){
            df[,i] = as.character(df[,i])
        }
    }
    if(qiimeVersion==TRUE){
        if(ncol(df)==0){
            meta = lapply(1:nrow(df),function(i){
                    ll = list(id=rownames(df)[i],metadata=NULL)
                    ll
                })
        } else {
        	meta = lapply(1:nrow(df),function(i){
            	    ll = list(id=rownames(df)[i], 
	                metadata=list("taxonomy" = paste(df[i,])))
                    NAvalues = grep("NA$",ll$metadata$taxonomy)
                    if(length(NAvalues)>0){
                        k = NAvalues[1]
                        ll$metadata$taxonomy = paste(df[i,1:(k-1)])
                    }
                	ll
            	})
        }
        return(meta)
    } else {
        if(ncol(df)==0){
            meta = lapply(1:nrow(df),function(i){
                    ll = list(id=rownames(df)[i],metadata=NULL)
                    ll
                })
        } else {
            meta  = lapply(1:nrow(df),function(i){
                    ll = list(id=rownames(df)[i],
                    metadata=lapply(1:ncol(df),
                        function(j){as.character(df[i,j])}))
                    names(ll$metadata) = colnames(df)
                    ll
                })            
        }
        return(meta)
    }
}
