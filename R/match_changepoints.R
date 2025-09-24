## GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' @title Match Detected and True Breakpoints
#' @description
#' Compares detected changepoints against the true changepoints using a tolerance
#' window. Each detected changepoint can be matched to at most one true changepoint,
#' and vice versa. This allows the computation of True Positives (TP),
#' False Positives (FP), and False Negatives (FN), as well as derived
#' performance metrics (precision, recall, F1-score).
#'
#' @param true_bkps An integer vector of true changepoint indices.
#' @param est_bkps An integer vector of detected changepoint indices.
#' @param w A non-negative integer specifying the tolerance window. A detected
#' changepoint within \code{w} indices of a true changepoint is considered a match.
#' Defaults to 5.
#'
#' @return A named list with the following elements:
#' \itemize{
#'   \item \code{TP}: Number of true positives (matched changepoints).
#'   \item \code{FP}: Number of false positives (unmatched detected changepoints).
#'   \item \code{FN}: Number of false negatives (unmatched true changepoints).
#'   \item \code{precision}: Precision score = TP / (TP + FP).
#'   \item \code{recall}: Recall score = TP / (TP + FN).
#'   \item \code{F1}: F1-score, harmonic mean of precision and recall.
#' }
#'
#' @examples
#' true_bkps <- c(100, 300, 400, 450, 550, 700, 750, 950, 1000)
#' est_bkps  <- c(98, 302, 398, 444, 555, 696, 742, 959, 997)
#'
#' # With tolerance window = 5, all estimated breakpoints match
#' match_breakpoints(true_bkps, est_bkps, w = 5)
#'
#' # With a stricter tolerance, some become FP/FN
#' match_breakpoints(true_bkps, est_bkps, w = 2)
#'
#' @export
match_breakpoints <- function(true_bkps, est_bkps, w = 5)
{
  true_bkps <- sort(true_bkps)
  est_bkps  <- sort(est_bkps)
  
  matched_true <- rep(FALSE, length(true_bkps))
  matched_est  <- rep(FALSE, length(est_bkps))
  
  for (i in seq_along(est_bkps))
  {
    diffs <- abs(true_bkps - est_bkps[i])
    candidates <- which(diffs <= w & !matched_true)
    if (length(candidates) > 0)
    {
      j <- candidates[which.min(diffs[candidates])]
      matched_true[j] <- TRUE
      matched_est[i]  <- TRUE
    }
  }
  
  TP <- sum(matched_true)
  FP <- sum(!matched_est)
  FN <- sum(!matched_true)
  
  precision <- ifelse(TP + FP > 0, TP / (TP + FP), NA)
  recall    <- ifelse(TP + FN > 0, TP / (TP + FN), NA)
  f1        <- ifelse(!is.na(precision) & !is.na(recall) & (precision+recall) > 0,
                      2 * precision * recall / (precision + recall), NA)
  
  list(TP = TP, FP = FP, FN = FN,
       precision = precision, recall = recall, F1 = f1)
}
