################################################################################
# metagenomeSeq plot functions unit tests
################################################################################
context("Testing norm factor calculation")
library("metagenomeSeq"); library("testthat")


test_that("`calcNormFactors` function provides expected values", {
  # uses the lung data and pre-calculated normalization factors 
  # for various values of p
  data(lungData)
  point25 = c(29,2475,2198,836,722,1820,79,1171,1985,710,145,742,848,89,1981)
  point = c(43,2475,2198,836,722,1820,119,1171,1985,710,145,742,848,89,1981)
  point100=as.numeric(unlist(libSize(lungData[,1:15])))
  expect_equal(as.numeric(unlist(calcNormFactors(lungData[,1:15]))),point)
  expect_equal(as.numeric(unlist(calcNormFactors(lungData[,1:15],p=.25))),point25)
  expect_equal(as.numeric(unlist(calcNormFactors(lungData[,1:15],p=1))),point100)
})

test_that("`cumNorm` returns the same object as defined in the package", {
  data(lungData); data(mouseData)
  expect_equal(cumNorm(mouseData,p=.5), mouseData)
  expect_equal(cumNorm(lungData), lungData)
})

test_that("`cumNormStat` returns the correct value", {
  data(lungData); data(mouseData);
  expect_equal(as.numeric(cumNormStat(lungData)),0.7014946)
  expect_equal(as.numeric(cumNormStat(mouseData)),0.5)
})

test_that("`cumNormStatFast` returns the correct value", {
  data(lungData); data(mouseData);
  expect_equal(as.numeric(cumNormStatFast(lungData)),0.7014946)
  expect_equal(as.numeric(cumNormStatFast(mouseData)),0.5)
})

