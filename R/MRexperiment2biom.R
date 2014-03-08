#' MRexperiment to biom objects
#' 
#' Wrapper to convert MRexperiment objects to biom objects.
#' 
#' @param obj The MRexperiment object.
#' @param id Optional id for the biom matrix.
#' @return A biom object.
#' @seealso \code{\link{load_meta}} \code{\link{load_phenoData}} \code{\link{newMRexperiment}} \code{\link{load_biom}} \code{\link{biom2MRexperiment}}
MRexperiment2biom <- function(obj,id=NULL){
	library(biom)
	id = id
	format = "Biological Observation Matrix 1.0.0-dev"
	format_url = "http://biom-format.org/documentation/format_versions/biom-1.0.html"
	type = "OTU table"
	generated_by = sprintf("metagenomeSeq %s",packageVersion("metagenomeSeq"))
	date = as.character(Sys.time())
	matrix_type = "dense"
	matrix_element_type = "int"

	data = MRcounts(obj)
	shape = dim(data)

	fdatanames = colnames(fData(obj))
	rows = lapply(1:nrow(data),function(i){
			ll = list(
				id=rownames(data)[i],
				metadata=lapply(1:ncol(fData(obj)),function(j){
					as.character(fData(obj)[i,j])}))
			names(ll$metadata) = fdatanames
			ll
			})

	sdatanames = colnames(pData(obj))
	columns  = lapply(1:ncol(data),function(i){
			ll = list(
				id=colnames(data)[i],
				metadata=lapply(1:ncol(pData(obj)),function(j){
					as.character(pData(obj)[i,j])}))
			names(ll$metadata) = sdatanames
			ll
			})
	data = as.list(as.data.frame(t(data)))
	names(data) <- NULL
	
	biomlist = list(id=id,format=format,format_url=format_url,type=type,generated_by=generated_by,
					date=date,matrix_type=matrix_type,matrix_element_type=matrix_element_type,shape=shape,
					rows=rows,columns=columns,data=data)
	biom(biomlist)
}
