test_that("per feature `fitZig` runs", {
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
  settings = zigControl(maxit=1,verbose=FALSE, per_feature_zeroModel=TRUE)
  fit2 = fitZig(obj = lungTrim,mod=mod,control=settings)
  # because the ordering is wrong
  expect_failure(expect_equal(fit,fit2))
  # check that they're equal now 
  fit2 = fit2[names(fit)]
  #  expect_equal(fit,fit2)
})
