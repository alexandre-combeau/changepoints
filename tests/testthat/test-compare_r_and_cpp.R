library(testthat)
library(changepoints)

test_that("smallest_valid_partitioning and smallest_valid_partitioning_rcpp produce similar results with valid_FOCUS", {
  set.seed(42)
  y <- c(rnorm(100), rnorm(100, 0.5))
  gamma <- 10

  # Use valid_FOCUS as the test function
  svp_r <- smallest_valid_partitioning(y, gamma = gamma, test = valid_FOCUS)
  svp_rcpp <- smallest_valid_partitioning_rcpp(y, gamma = gamma, test = valid_FOCUS)

  # Their changepoints should be similar
  expect_equal(svp_r$changepoints, svp_rcpp$changepoints, tolerance = 0)

  # Their costQ vectors should be nearly identical
  expect_equal(svp_r$costQ, svp_rcpp$costQ, tolerance = 1e-6)
})