##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' @title Evaluate Changepoint Detection Accuracy
#' @description
#' Computes precision, recall, and F1-score for a set of estimated changepoints
#' compared to ground truth, within a user-defined tolerance.
#'
#' @param algo A character string giving the name of the algorithm being evaluated.
#'   Defaults to `"MyAlgo"`.
#' @param cp_true An integer vector (or list) of true changepoint indices.
#' @param cp_est An integer vector (or list) of estimated changepoint indices.
#' @param tol A non-negative integer giving the tolerance window for matching
#'   estimated changepoints to true changepoints. Defaults to `5`.
#'
#' @details
#' The function matches estimated changepoints to true ones using a greedy
#' nearest-neighbor approach within the specified tolerance:
#' \itemize{
#'   \item True Positives (TP): correctly identified changepoints.
#'   \item False Positives (FP): spurious estimated changepoints.
#'   \item False Negatives (FN): missed true changepoints.
#' }
#'
#' Precision, recall, and F1 are then computed as:
#' \deqn{Precision = TP / (TP + FP)}
#' \deqn{Recall = TP / (TP + FN)}
#' \deqn{F1 = 2 * Precision * Recall / (Precision + Recall)}
#'
#' If no changepoints are present in both `cp_true` and `cp_est`, all metrics are returned as `NA`.
#' If only false alarms occur (`cp_true` empty), precision is `0` and recall is `NA`.
#' If only misses occur (`cp_est` empty), recall is `0` and precision is `NA`.
#'
#' @return A data frame with one row and the following columns:
#' \itemize{
#'   \item `Algo`: the algorithm name,
#'   \item `Precision`: precision score,
#'   \item `Recall`: recall score,
#'   \item `F1`: F1-score.
#' }
#'
#' @examples
#' # Example with true changepoints at positions 50 and 100
#' cp_true <- c(50, 100)
#' 
#' # Estimated changepoints, one correct within tolerance, one false alarm
#' cp_est <- c(52, 150)
#' 
#' cp_metrics("MyAlgo", cp_true, cp_est, tol = 5)
#'
#' # Case with perfect detection
#' cp_metrics("PerfectAlgo", cp_true, cp_true, tol = 5)
#'
#' # Case with no estimates
#' cp_metrics("EmptyAlgo", cp_true, integer(0), tol = 5)
#'
#' @export
cp_metrics <- function(algo = "MyAlgo", cp_true, cp_est, tol = 5)
{
  # Force plain integer vectors (handles list-columns safely)
  cp_true <- as.integer(unlist(cp_true))
  cp_est  <- as.integer(unlist(cp_est))
  
  # Edge cases
  if (length(cp_true) == 0 && length(cp_est) == 0)
  {
    data.frame(
      Algo = algo,
      Precision = NA_real_,
      Recall = NA_real_,
      F1 = NA_real_
    )
  }
  
  if (length(cp_true) == 0) # only false alarms
  {
    precision <- 0
    recall    <- NA_real_
    
    data.frame(
      Algo = algo,
      Precision = precision,
      Recall = recall,
      F1 = NA_real_
    )
  }
  if (length(cp_est) == 0) # all misses
  {
    precision <- NA_real_
    recall    <- 0
    data.frame(
      Algo = algo,
      Precision = precision,
      Recall = recall,
      F1 = NA_real_
    )
  }
  
  # Pairwise absolute distances
  D <- abs(outer(cp_true, cp_est, "-"))
  
  # Greedy matching: just pick the closest unmatched estimate within tolerance
  matched_est <- rep(FALSE, length(cp_est))
  TP <- 0
  for (i in seq_along(cp_true))
  {
    j <- which.min(D[i, ])
    if (length(j) && D[i, j] <= tol && !matched_est[j])
    {
      TP <- TP + 1
      matched_est[j] <- TRUE
    }
  }
  
  FP <- sum(!matched_est)
  FN <- length(cp_true) - TP
  
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall) else NA_real_
  
  data.frame(
    Algo = algo,
    Precision = precision,
    Recall = recall,
    F1 = f1
  )
}
