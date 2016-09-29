#' Load objects organized in the Biom format.
#' 
#' Wrapper to load Biom formatted object. 
#' 
#' @param file The biom object filepath.
#' @return A MRexperiment object.
#' @seealso \code{\link{load_meta}} \code{\link{load_phenoData}} \code{\link{newMRexperiment}} \code{\link{biom2MRexperiment}}
#' @examples
#' 
#' #library(biomformat)
#' #rich_dense_file = system.file("extdata", "rich_dense_otu_table.biom", package = "biomformat")
#' #x = loadBiom(rich_dense_file)
#' #x
loadBiom <- function(file){
	requireNamespace("biomformat")
	x = biomformat::read_biom(file);
	mrobj = biom2MRexperiment(x);
	return(mrobj);
}
