#' Shrink coefficients of a fitZig result
#' 
#' @param obj result of call to `fitZig` or `.do_fitZig` functions
#' @param coef coefficient to shrink
#' 
#' @result list with same structure as input object with coefficient updated
#' 
shrinkZig <- function(obj, coef=2) {
  eb <- obj$eb
  y <- obj$y
  
  if (is.null(y)) {
    if (!is.null(obj$counts)) {
      y <- obj$counts
    } else {
      stop("Need y matrix to shrink coefficient")
    }
  } 

  # calcuate variance of coefficient tau
  betas <- eb$coefficients[,coef]
  tau <- var(betas, na.rm=TRUE) * nrow(y)

  # get residuals with reduced model
  mat1 <- eb$design
  coefs <- eb$coefficients
  
  if (ncol(coefs) > 2) {
    coef_ind <- c(1,coef)
    mat0 <- mat1[,-coef_ind,drop=FALSE]
    mat1 <- mat1[,coef_ind,drop=FALSE]
    
    coefs0 <- coefs[, -coef_ind]
    ybar <- coefs0 %*% t(mat0)
    y <- y - ybar
  }
  
  # fit ridge regression for coef
  shrunk_coef <- sapply(seq_len(nrow(y)), function(i) {
    zi <- 1-obj$z[i,]
    lambda <- eb$sigma[i]^2 / (sum(zi) * tau)
    ridge_res <- glmnet(y=y[i,], x=mat1, 
                        lambda=lambda, alpha=0, 
                        weights=zi, 
                        intercept=TRUE,
                        standardize=FALSE)
    b0 <- coefficients(ridge_res)[3,1]
    b0
  })
  
  # update coef in object, return obj
  eb$coefficients[,coef] <- shrunk_coef
  obj$eb <- eb
  return(obj)
}