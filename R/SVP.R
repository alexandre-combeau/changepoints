##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' Smallest Valid Partitioning with Validation and Pruning
#' 
#' @title Smallest Valid Partitioning with Validation and Pruning
#' @description
#' This function implements a dynamic programming approach to segment a univariate signal 
#' into the smallest number of valid segments, according to a user-defined validation function. 
#' Each segment must pass a validity test (e.g., based on variance, range, etc.). 
#' The algorithm minimizes a quadratic cost subject to this constraint.
#'
#' It uses cumulative sums for efficient cost calculation, and optionally applies pruning 
#' to discard changepoint candidates that would lead to invalid segments.
#'
#' @param data A numeric vector representing the univariate signal to be segmented.
#' @param gamma A numeric value used as a threshold in the validation function and as a penalty for each segment.
#' @param test A function of the form `function(data, gamma)` returning TRUE if the segment is valid. Default is `valid_OP`.
#' @param all_full_validity Logical. If TRUE (default), the algorithm applies *segment-wise validation*:
#' at each time step, it tests whether the candidate segment \code{data[(s+1):t]} is valid using the 
#' user-defined function \code{test}. If the segment fails the test, the candidate \code{s} is removed 
#' (pruned) from the set of possible changepoints. This accelerates computation by avoiding invalid 
#' segment extensions. If FALSE, the algorithm skips this validation and considers all candidate 
#' segments without checking their validity (which can be faster but may return invalid segments).
#' 
#' @return A list with the following components :
#' \describe{
#'   \item{changepoints}{Integer vector indicating the ending index of each segment (i.e., positions of changepoints).}
#'   \item{nb}{Integer vector of length \code{length(data)}. At each position \code{t}, it records the number of candidates tested.}
#'   \item{costQ}{Numeric vector of length \code{length(data)}. Quadratic cost value at each time step.}
#'   \item{R}{A matrix of dimension \code{length(data) x 3} containing, for each time step: 
#'     \code{Q} (cumulative cost), \code{K} (number of segments), and \code{s} (previous changepoint).}
#' }
#' 
#' @examples
#' data <- c(rnorm(100, 0, 1), rnorm(100, 5, 1))
#' SVP <- smallest_valid_partitioning(data, gamma = 10, test = valid_OP)
#' plot_segmentation(data, SVP$changepoints, title = "SVP_OP Segmentation")
#'
#' @export
smallest_valid_partitioning <- function(data, gamma, test = valid_OP, all_full_validity = TRUE)
{
  n <- length(data)
  
  # Initialisation de la matrice (Q, K, s)
  R <- matrix(Inf, nrow = n+1, ncol = 3)
  colnames(R) <- c("Q", "K", "s")
  
  # Condition initiale t=0
  R[1, ] <- c(0, 0, 0)
  
  nb <- integer(n)    # nb de candidats examinés à chaque t
  costQ <- numeric(n) # Q_t final pour chaque t
  
  # Cumulative sum for optimized calculations
  cs_y  <- c(0, cumsum(data))
  cs_y2 <- c(0, cumsum(data^2))
  
  INDEX <- 0
  
  for (t in 1:n)
  {
    best_Q <- Inf         # best cost for the current t
    best_K <- Inf         # best segment number
    best_s <- NA_integer_ # last changepoint

    for (s in INDEX)
    {
      if (all_full_validity) # tester si tous les sous-segments passent le test de f<gamma
      {
        if (!test(data[(s+1):t], gamma))
        {
          INDEX <- setdiff(INDEX, s) # we remove the invalid candidate
          next # skipping the iteration
        }
      }
      
      segment_cost <- (cs_y2[t+1] - cs_y2[s+1]) - (cs_y[t+1] - cs_y[s+1])^2 / (t-s)
      candidate_Q <- R[s+1, "Q"] + segment_cost
      candidate_K <- R[s+1, "K"] + 1
      
      if (candidate_K < best_K || (candidate_K == best_K && candidate_Q < best_Q))
      {
        best_Q <- candidate_Q
        best_K <- candidate_K
        costQ[t] <- best_Q
        best_s <- s
      }
    }
    
    for (s in INDEX)
    {
      if (s != t)
      {
        segment_cost <- (cs_y2[t+1]-cs_y2[s+1]) - (cs_y[t+1]-cs_y[s+1])^2 / (t-s)
        inequality_test <- R[s + 1, "Q"] + segment_cost > best_Q
        equality_test   <- R[s + 1, "K"] == best_K
        if (inequality_test && equality_test)
        {
          INDEX <- setdiff(INDEX, s)
        }
      }
    }
    
    nb[t] <- length(INDEX)
    INDEX <- c(INDEX, t)
    
    # we keep the best candidate for R_t
    R[t+1, ] <- c(best_Q, best_K, best_s)
  }
  
  # Reconstruction
  changepoints <- integer(0)
  t <- n
  while (t > 0)
  {
    changepoints <- c(t, changepoints) # last index of the current segment
    t <- R[t + 1, "s"]                 # back to the previous s for the next loop
  }
  
  list(changepoints = unname(changepoints),
       lastIndexSet = rev(INDEX),
       nb           = nb,
       costQ        = costQ,
       R            = R[-1,])
}
