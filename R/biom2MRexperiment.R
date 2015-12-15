#' Biom to MRexperiment objects
#' 
#' Wrapper to convert biom files to MRexperiment objects.
#' 
#' @param obj The biom object file.
#' @return A MRexperiment object.
#' @seealso \code{\link{load_meta}} \code{\link{load_phenoData}} \code{\link{newMRexperiment}} \code{\link{load_biom}}
#' @examples
#' 
#' #library(biomformat)
#' #rich_dense_file = system.file("extdata", "rich_dense_otu_table.biom", package = "biomformat")
#' #x = read_biom(rich_dense_file)
#' #biom2MRexperiment(x)
biom2MRexperiment <- function(obj){
	requireNamespace("biomformat")
	mat = as(biomformat::biom_data(obj),"matrix")

	if(! is.null(biomformat::observation_metadata(obj))){
		len = max(sapply(biomformat::observation_metadata(obj),length))
		taxa = as.matrix(sapply(biomformat::observation_metadata(obj),function(i){ i[1:len]}))
		
		if(dim(taxa)[1]!=dim(mat)[1]){
			taxa = t(taxa)
		}
		rownames(taxa) = rownames(mat)
		colnames(taxa) = colnames(biomformat::observation_metadata(obj))
		taxa = as(data.frame(taxa),"AnnotatedDataFrame")
	} else{
		taxa = NULL
	}

	if(! is.null(biomformat::sample_metadata(obj))) {
		pd = as(biomformat::sample_metadata(obj),"AnnotatedDataFrame")
	} else{
		pd = NULL
	}
	
	mrobj = newMRexperiment(counts = mat, phenoData = pd, featureData = taxa)
	return(mrobj)
}
