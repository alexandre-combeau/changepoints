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
  # Their lastIndexSet should be similar
  expect_equal(svp_r$lastIndexSet, svp_rcpp$lastIndexSet, tolerance = 0)
  # Their nb of candidates examinated should be similar
  expect_equal(svp_r$nb, svp_rcpp$nb, tolerance = 0)

  # Their R matrix should be similar
  expect_equal(svp_r$R[,1], svp_rcpp$R[,1], tolerance = 1e-6)
  expect_equal(svp_r$R[,2], svp_rcpp$R[,2], tolerance = 0)
  expect_equal(svp_r$R[,3], svp_rcpp$R[,3], tolerance = 0)
})


test_that("smallest_valid_partitioning_rcpp_seq and smallest_valid_partitioning_rcpp produce similar results for gaussian_mean", {
  set.seed(123)
  y <- c(rnorm(100), rnorm(100, 0.5))
  gamma <- 1
  
  # Classic Rcpp implementation (with valid_FOCUS)
  svp_rcpp <- smallest_valid_partitioning_rcpp(y, gamma, test = valid_FOCUS,
                                               prune_if_unvalid = T,
                                               prune_if_PELT = T)
  
  # Sequential implementation (with gaussian_mean)
  svp_seq <- smallest_valid_partitioning_rcpp_seq_new(y, gamma, test = "gaussian_mean",F,T)
  
  # Their changepoints should be similar
  expect_equal(svp_rcpp$changepoints, svp_seq$changepoints, tolerance = 0)
  # Their lastIndexSet should be similar
  expect_equal(svp_rcpp$lastIndexSet, svp_seq$lastIndexSet, tolerance = 0)
  # Their nb of candidates examinated should be similar
  expect_equal(svp_rcpp$nb, svp_seq$nb, tolerance = 0)

  # Their R matrix should be similar
  expect_equal(svp_rcpp$R[,1], svp_seq$R[,1], tolerance = 1e-6)
  expect_equal(svp_rcpp$R[,2], svp_seq$R[,2], tolerance = 0)
  expect_equal(svp_rcpp$R[,3], svp_seq$R[,3], tolerance = 0)
})


test_that("smallest_valid_partitioning_rcpp_seq and smallest_valid_partitioning_rcpp produce same results (no pruning)", {
  set.seed(123)
  y <- c(rnorm(500), rnorm(10, 0.5))
  gamma <- 100

  # Classic Rcpp implementation (with valid_FOCUS)
  svp_rcpp <- smallest_valid_partitioning_rcpp(y, gamma, test = valid_FOCUS)

  # Sequential implementation (with gaussian_mean)
  svp_seq <- smallest_valid_partitioning_rcpp_seq_old(y, gamma, test = "gaussian_mean",
                                                  all_full_validity = TRUE)

  # Their changepoints should be similar
  expect_equal(svp_rcpp$changepoints, svp_seq$changepoints, tolerance = 0)
  # Their nb of candidates examinated should be similar
  expect_equal(svp_rcpp$nb, svp_seq$nb, tolerance = 0)

  # Their R matrix should be similar
  expect_equal(svp_rcpp$R[,1], svp_seq$R[,1], tolerance = 1e-6)
  expect_equal(svp_rcpp$R[,2], svp_seq$R[,2], tolerance = 0)
  expect_equal(svp_rcpp$R[,3], svp_seq$R[,3], tolerance = 0)
})



##### TEST BETWEEN SEQ_OLD AND SEQ_NEW #####
test_that("smallest_valid_partitioning_rcpp_seq_old and smallest_valid_partitioning_rcpp_seq_new should produce same results", {
  set.seed(123)
  y <- c(rnorm(100), rnorm(100, 0.5))
  gamma <- 1
  
  # Classic Rcpp implementation (with valid_FOCUS)
  svp_seq_old <- smallest_valid_partitioning_rcpp_seq_old(y, gamma, test = "gaussian_mean",
                                                          all_full_validity = TRUE)
  
  # Sequential implementation (with gaussian_mean)
  svp_seq_new <- smallest_valid_partitioning_rcpp_seq_new(y, gamma, test = "gaussian_mean",
                                                          prune_if_unvalid = FALSE,
                                                          prune_if_PELT = TRUE)
  
  # Their changepoints should be similar
  expect_equal(svp_seq_old$changepoints, svp_seq_new$changepoints, tolerance = 0)
  # Their nb of candidates examinated should be similar
  expect_equal(svp_seq_old$nb, svp_seq_new$nb, tolerance = 0)
  
  # Their R matrix should be similar
  expect_equal(svp_seq_old$R[,1], svp_seq_new$R[,1], tolerance = 1e-6)
  expect_equal(svp_seq_old$R[,2], svp_seq_new$R[,2], tolerance = 0)
  expect_equal(svp_seq_old$R[,3], svp_seq_new$R[,3], tolerance = 0)
})
