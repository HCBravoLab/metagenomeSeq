#' Load objects organized in the Biome format.
#' 
#' Wrapper to load Biome formatted object. 
#' 
#' @param file The biome object filepath.
#' @return A MRexperiment object.
#' @seealso \code{\link{load_meta}} \code{\link{load_phenoData}} \code{\link{newMRexperiment}} \code{\link{biom2MRexperiment}}
#' @examples
#' 
#' #library(biom)
#' #rich_dense_file = system.file("extdata", "rich_dense_otu_table.biom", package = "biom")
#' #x = load_biome(rich_dense_file)
#' #x
load_biom <- function(file){
	library(biom)
	x = biom::read_biom(file);
	mrobj = biom2MRexperiment(x);
	return(mrobj);
}
