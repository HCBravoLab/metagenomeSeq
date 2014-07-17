#' Correlation of each row of a matrix or MRexperiment object
#'
#' Calculates the (pairwise) correlation statistics and associated p-values of a matrix
#' or the correlation of each row with a vector.
#'
#' @param obj A MRexperiment object or count matrix.
#' @param y Vector of length ncol(obj) to compare to.
#' @param method One of 'pearson','spearman', or 'kendall'.
#' @param alternative Indicates the alternative hypothesis and must be one of 'two.sided', 'greater' (positive) or 'less'(negative). You can specify just the initial letter.
#' @param norm Whether to aggregate normalized counts or not - if MRexperiment object.
#' @param log Whether or not to log2 transform the counts - if MRexperiment object.
#' @param parallel Parallelize the correlation testing?
#' @param cores Number of cores to make use of if parallel == TRUE.
#' @param override If the number of rows to test is over a thousand the test will not commence (unless override==TRUE).
#' @param ... Extra settings for mclapply.
#' @return A matrix of size choose(number of rows, 2) by 2. The first column corresponds to the correlation value. The second column the p-value.
#' @seealso \code{\link{correctIndices}}
#' @aliases corTest
#' @export
#' @examples
#'
#' # Pairwise correlation of raw counts
#' data(mouseData)
#' cors = correlationTest(mouseData[1:10,],libSize(mouseData),norm=FALSE,log=FALSE)
#' head(cors)
#' 
#' cormat = as.matrix(dist(mat)) # Creating a matrix
#' cormat[cormat>0] = 0 # Creating an empty matrix
#' cormat[upper.tri(cormat)][ind] = cors[,1] 
#' table(cormat[1,-1] - cors[1:5,1])
#'
#' # Correlation of raw counts with a vector (library size in this case)
#' data(mouseData)
#' cors = correlationTest(mouseData[1:10,],libSize(mouseData),norm=FALSE,log=FALSE)
#' head(cors)
#'
correlationTest <- function(obj,y=NULL,method="pearson",alternative="two.sided",norm=TRUE,log=TRUE,parallel=FALSE,cores=2,override=FALSE,...){
	if(class(obj)=="MRexperiment"){
		mat = MRcounts(obj,norm=norm,log=log)
	} else if(class(obj) == "matrix") {
		mat = obj
	} else {
		stop("Object needs to be either a MRexperiment object or matrix")
	}
	nr = nrow(mat)
	if(nr > 1000){
		if(override){
			show("Good luck! This might take some time.")
		} else {
			stop("Many features being considered - to proceed set override to TRUE")
		}
	}

	if(is.null(rownames(mat))){
		nm = as.character(1:nr)
	} else { 
		nm = rownames(mat)
	}

	if(parallel){
		if(require(parallel)){
			if(is.null(y)){
				corrAndP = mclapply(1:(nr-1),function(i){
					vals =(i+1):nr
					cp = array(NA,dim=c(length(vals),2))
					rownames(cp) = paste(nm[i],nm[(i+1):nr],sep="-")
					colnames(cp) = c("correlation","pvalue")
					for(j in (i+1):nr){	
						x = as.numeric(mat[i,])
						y = as.numeric(mat[j,])
						res = cor.test(x,y,method=method,
							alternative=alternative)
						cp[j-i,1] = res$estimate
						cp[j-i,2] = res$p.value
					}
					cp
				},mc.cores=cores,...)
			} else {
				corrAndP = mclapply(1:nr,function(i){
					res = cor.test(mat[i,],y,method=method,
						alternative=alternative)
					cbind(res$estimate,res$p.value)
				},mc.cores=cores,...)
			}
		} else {
			parallel=FALSE
		}
	}
	if(parallel==FALSE){
		if(is.null(y)){
			corrAndP=lapply(1:(nr-1),function(i){			
				vals =(i+1):nr
				cp = array(NA,dim=c(length(vals),2))
				rownames(cp) = paste(nm[i],nm[(i+1):nr],sep="-")
				colnames(cp) = c("correlation","pvalue")
				for(j in (i+1):nr){	
					x = as.numeric(mat[i,])
					y = as.numeric(mat[j,])
					res = cor.test(x,y,method=method,
						alternative=alternative)
					cp[j-i,1] = res$estimate
					cp[j-i,2] = res$p.value
				}
				cp
			})
		} else {
			corrAndP=lapply(1:nr,function(i){
					res = cor.test(mat[i,],y,method=method,
						alternative=alternative)
					cbind(res$estimate,res$p.value)
				})
		}
	}
	# browser()
	correlation = unlist(sapply(corrAndP,function(i){i[,1]}))
	p  = unlist(sapply(corrAndP,function(i){i[,2]}))
	results = cbind(correlation,p)
	if(!is.null(y)) rownames(results) = rownames(obj)
	
	return(results)
}
#' Calculate the correct indices for the output of correlationTest
#'
#' Consider the upper triangular portion of a matrix of size nxn. Results from the \code{correlationTest} are output
#' as the combination of two vectors, correlation statistic and p-values. The order of the output is 1vs2, 1vs3, 1vs4, etc.
#' The correctIndices returns the correct indices to fill a correlation matrix or correlation-pvalue matrix.
#'
#' @param n The number of features compared by correlationTest (nrow(mat)).
#' @return A vector of the indices for an upper triangular matrix.
#' @seealso \code{\link{correlationTest}}
#' @export
#' @examples
#' 
#' data(mouseData)
#' mat = MRcounts(mouseData)[55:60,]
#' cors = correlationTest(mat)
#' ind  = correctIndices(nrow(mat))
#' 
#' cormat = as.matrix(dist(mat))
#' cormat[cormat>0] = 0
#' cormat[upper.tri(cormat)][ind] = cors[,1]
#' table(cormat[1,-1] - cors[1:5,1])
#'
correctIndices <- function(n){
	if(n==1){
		return(1)
	}
	if(n==2){
		return(c(1,2))
	}
	seq1 <- cumsum(1:(n-1)) - c(0,1:(n-2))
	seq2 <- sapply(1:(n-2),function(i) {
		seq1[-c(1:i)]+1*i  
	})
	seq <- c(seq1,unlist(seq2))
	return(seq)
}