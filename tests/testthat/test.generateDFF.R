context("Distinguishing Features")
library(phydose)


test_that("generateDFF returns the correct data structure", {
  expect_equal(is.list(generateDistFeat(phydata$trees)), TRUE)
  expect_equal(is.list(generateDistFeat(phydata$trees)[[1]]), TRUE)
  expect_equal(is.list(generateDistFeat(phydata$trees)[[1]][[1]]), FALSE)
  
  
})

