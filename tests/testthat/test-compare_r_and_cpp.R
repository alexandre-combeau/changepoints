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


test_that("smallest_valid_partitioning_rcpp_seq and smallest_valid_partitioning_rcpp produce similar results for gaussian_mean", {
  set.seed(42)
  y <- c(rnorm(100), rnorm(100, 0.5))
  gamma <- 1

  # Classic Rcpp implementation (with valid_FOCUS)
  svp_rcpp <- smallest_valid_partitioning_rcpp(y, gamma = gamma, test = valid_FOCUS)

  # Sequential implementation (with gaussian_mean)
  svp_seq <- smallest_valid_partitioning_rcpp_seq(y, gamma = gamma, test = "gaussian_mean", all_full_validity = TRUE)

  # Compare changepoints
  expect_equal(svp_seq$changepoints, svp_rcpp$changepoints, tolerance = 0)

  # Compare costQ vectors
  expect_equal(svp_seq$costQ, svp_rcpp$costQ, tolerance = 1e-6)
})


test_that("smallest_valid_partitioning_rcpp_seq and smallest_valid_partitioning_rcpp produce same results (no pruning)", {
  set.seed(123)
  y <- c(rnorm(500), rnorm(10, 0.5))
  gamma <- 100

  # Classic Rcpp implementation (with valid_FOCUS)
  svp_rcpp <- smallest_valid_partitioning_rcpp(y, gamma = gamma, test = valid_FOCUS)

  # Sequential implementation (with gaussian_mean)
  svp_seq <- smallest_valid_partitioning_rcpp_seq(y, gamma = gamma, test = "gaussian_mean", all_full_validity = TRUE)

  # Compare changepoints
  expect_equal(svp_seq$changepoints, svp_rcpp$changepoints, tolerance = 0)

  # Compare costQ vectors
  expect_equal(svp_seq$costQ, svp_rcpp$costQ, tolerance = 1e-6)
})