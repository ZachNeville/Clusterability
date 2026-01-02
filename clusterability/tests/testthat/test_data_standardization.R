context("Data Standardization")
library(clusterability)

test_that("standardize_data", {
  # Setup
  constcol <- matrix(c(1, 1, 1, 4, 6, 10, 9, 5, 3), nrow = 3)

  set.seed(1234)
  data2 <- cbind(rnorm(100, 7, 2), rnorm(100, 4, 3))

  noneresult <- standardize_data(data2, "NONE")
  stdresult <- standardize_data(data2, "STD")
  meanresult <- standardize_data(data2, "MEAN")
  medianresult <- standardize_data(data2, "MEDIAN")

  # Tests - Error handling
  expect_warning(standardize_data(constcol, "STD"), info = "Warn when constant variable causes NaN's")

  # Tests - is the output actually correct?
  expect_equal(colMeans(noneresult), colMeans(data2), tolerance = 1e-14)
  expect_equal(colMeans(stdresult), c(0, 0), tolerance = 1e-14)
  expect_equal(colMeans(meanresult), c(0, 0), tolerance = 1e-14)
  expect_equal(apply(medianresult, 2, median), c(0, 0), tolerance = 1e-14)

  expect_equal(apply(noneresult, 2, sd), apply(data2, 2, sd), tolerance = 1e-14)
  expect_equal(apply(stdresult, 2, sd), c(1, 1), tolerance = 1e-14)
  expect_equal(apply(meanresult, 2, sd), apply(data2, 2, sd), tolerance = 1e-14)
  expect_equal(apply(medianresult, 2, sd), apply(data2, 2, sd), tolerance = 1e-14)
})
