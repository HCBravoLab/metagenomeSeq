## unit test for wrenchNorm

context("Test that wrenchNorm functions properly")
library("metagenomeSeq"); library("testthat");

test_that( "norm factors generated are correct",{
  data("lungData"); data("mouseData");
  mouseNF <- c(0.3364660,0.7051424,1.3295084,0.8530978,0.7545386,2.1273695,1.2158941,1.9025748,0.5382427,0.5841864)
  lungNF <- c(0.006551719,12.267861013,10.106967942,2.447679975,1.266012939,5.701245412,0.049474404,2.863477065,6.821474324,1.261155349)
  lungData <- lungData[, -which(is.na(pData(lungData)$SmokingStatus))]
  lungData2 <- wrenchNorm(lungData, condition = lungData$SmokingStatus)
  mouseData2 <- wrenchNorm(mouseData, condition = mouseData$diet)
  expect_equal(as.numeric(normFactors(lungData2)[1:10]), lungNF)
  expect_equal(as.numeric(unlist(normFactors(mouseData2)[1:10])), mouseNF, tolerance = 1e-06)
})
