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
	mat = as(biom_data(obj),"matrix")
	taxa = as(observation_metadata(obj),"AnnotatedDataFrame")
	pd = as(sample_metadata(obj),"AnnotatedDataFrame")
	mrobj = newMRexperiment(counts = mat, phenoData = pd, featureData = taxa);
	return(mrobj);
}
