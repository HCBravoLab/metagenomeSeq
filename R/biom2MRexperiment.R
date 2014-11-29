#' Biome to MRexperiment objects
#' 
#' Wrapper to convert biome files to MRexperiment objects.
#' 
#' @param obj The biome object file.
#' @return A MRexperiment object.
#' @seealso \code{\link{load_meta}} \code{\link{load_phenoData}} \code{\link{newMRexperiment}} \code{\link{load_biom}}
#' @examples
#' 
#' #library(biom)
#' #rich_dense_file = system.file("extdata", "rich_dense_otu_table.biom", package = "biom")
#' #x = read_biom(rich_dense_file)
#' #biom2MRexperiment(x)
biom2MRexperiment <- function(obj){
	library(biom)
	mat = as(biom::biom_data(obj),"matrix")

	if(! is.null(biom::observation_metadata(obj))){
		len = max(sapply(biom::observation_metadata(obj),length))
		cnames = names(biom::observation_metadata(obj)[[which.max(sapply(biom::observation_metadata(obj),length))]])
		taxa = as.matrix(sapply(biom::observation_metadata(obj),function(i){ i[1:len]}))
		
		if(dim(taxa)[1]!=dim(mat)[1]){
			taxa = t(taxa)
		}
		rownames(taxa) = rownames(mat)
		colnames(taxa) = cnames
		taxa = as(data.frame(taxa),"AnnotatedDataFrame")
	} else{
		taxa = NULL
	}

	if(! is.null(biom::sample_metadata(obj))) {
		pd = as(biom::sample_metadata(obj),"AnnotatedDataFrame")
	} else{
		pd = NULL
	}
	
	mrobj = newMRexperiment(counts = mat, phenoData = pd, featureData = taxa)
	return(mrobj)
}
