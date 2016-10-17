################################################################################
# metagenomeSeq plot functions unit tests
################################################################################
library("metagenomeSeq"); library("testthat"); 

test_that("`fitZig` function provides expected values prior to split", {
  # uses the lung data and pre-calculated fitZig result from 
  # prior to this separation
  data(lungData)
  path = system.file("extdata", package = "metagenomeSeq")
  fit  = readRDS(file.path(path,"lungfit.rds"))
  
  # run the same fit
  k = grep("Extraction.Control",pData(lungData)$SampleType)
  lungTrim = lungData[,-k]
  k = which(rowSums(MRcounts(lungTrim)>0)<30)
  lungTrim = cumNorm(lungTrim)
  lungTrim = lungTrim[-k,]
  smokingStatus = pData(lungTrim)$SmokingStatus
  mod = model.matrix(~smokingStatus)
  settings = zigControl(maxit=1,verbose=FALSE)
  fit2 = fitZig(obj = lungTrim,mod=mod,control=settings)
  # because the ordering is wrong
  expect_failure(expect_equal(fit,fit2))
  # check that they're equal now 
  fit2 = fit2[names(fit)]
  expect_equal(fit,fit2)
})

test_that("`fitZig` function treats a matrix the same", {
  # uses the lung data and pre-calculated fitZig result from 
  # prior to this separation
  data(lungData)
  path = system.file("extdata", package = "metagenomeSeq")
  fit  = readRDS(file.path(path,"lungfit.rds"))
  
  # run the same fit
  k = grep("Extraction.Control",pData(lungData)$SampleType)
  lungTrim = lungData[,-k]
  k = which(rowSums(MRcounts(lungTrim)>0)<30)
  lungTrim = cumNorm(lungTrim)
  lungTrim = lungTrim[-k,]
  smokingStatus = pData(lungTrim)$SmokingStatus
  scalingFactor = log2(normFactors(lungTrim)/1000 +1)
  mod = model.matrix(~smokingStatus)
  mod = cbind(mod,scalingFactor)
  settings = zigControl(maxit=1,verbose=FALSE)
  cnts = MRcounts(lungTrim)
  fit2 = fitZig(obj = lungTrim,mod=mod,control=settings,useCSSoffset=FALSE)
  fit2 = fit2[names(fit)]

  # expecting failure because of call
  expect_failure(expect_equal(fit,fit2))
  fit2$call = "123"
  fit$call = "123"
  # check that they're equal
  expect_equal(fit,fit2)
})
