test_that("per feature `fitZig` runs with shrinkage", {
  # uses the lung data and pre-calculated fitZig result from 
  # prior to this separation
  data(lungData)

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
  fit2 <- shrinkZig(fit2, coef=2)
  expect_is(fit2, "list")
})
