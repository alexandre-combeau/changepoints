##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' Validity Test Based on FOCuS
#' 
#' @title Validity Test Based on FOCuS
#' @description Checks whether a given segment is valid using the FOCuS method (i.e., has no changepoint).
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold passed to the FOCuS method.
#' @return A logical value indicating whether the segment is considered valid (TRUE if no changepoint detected).
#' @export
valid_focus <- function(y, gamma)
{
  res <- FOCuS(y, gamma)
  
  if(res$changepoint == -1)
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}

################################################

#' Validity Test Based on Variance
#' 
#' @title Validity Test Based on Variance
#' @description Checks if the variance of a segment is lower than or equal to a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the maximum allowed variance.
#' @return TRUE if the variance is less than or equal to gamma, or if the segment is very short.
#' @export
valid_var <- function(y, gamma)
{
  if(length(y) < 6)
  {
    return(TRUE)
  }
  
  sum((y - mean(y))^2) <= gamma
}

################################################

#' Validity Test Based on Range
#' 
#' @title Validity Test Based on Range
#' @description Checks whether the range (max - min) of the segment is below a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the maximum allowed range.
#' @return TRUE if the range is less than or equal to gamma.
#' @export
valid_Range <- function(y, gamma)
{
  max(y) - min(y) <= gamma
}

################################################

#' Validity Test Based on Trimmed Range
#' 
#' @title Validity Test Based on Trimmed Range
#' @description Applies a slack range test using trimmed minimum and maximum (ignores smallest and largest 3 values).
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the trimmed range.
#' @return TRUE if the trimmed range is less than or equal to gamma, or if the segment is too short.
#' @export
valid_Range_slack <- function(y, gamma)
{
  if(length(y) < 7)
  {
    return(TRUE)
  }
  sortedY <- sort(y)
  sortedY[length(y) - 3] - sortedY[3] <= gamma
}

################################################

#' Validity Test Based on Interquantile Range
#' 
#' @title Validity Test Based on Interquantile Range
#' @description Checks whether the interquantile range (95% - 5%) is below a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the interquantile range.
#' @return TRUE if the interquantile range is less than or equal to gamma.
#' @export
valid_Quantile <- function(y, gamma)
{
  quantile(y, probs = 0.95) - quantile(y, probs = 0.05) <= gamma
}

################################################
#' Optimal Partitioning Cost Test
#' 
#' @title Optimal Partitioning Cost Test
#' @description Tests whether the total cost of the segment is close enough to the best two-part segmentation cost.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for allowable difference between total cost and optimal bipartition.
#' @return TRUE if the segment cannot be split into two parts with significantly lower cost.
#' @export
valid_OP <- function(y, gamma)
{
  len <- length(y)
  if(len == 1){return(TRUE)}

  total <- sum((y - mean(y))^2)
  val <- Inf
  for(i in 2:len)
  {
    segment1 <- y[1:(i-1)]
    segment2 <- y[i:len]
    mean_seg1 <- mean(segment1)
    mean_seg2 <- mean(segment2)
    temp <- (sum((segment1 - mean_seg1)^2) + sum((segment2 - mean_seg2)^2))
    if(temp < val){val <- temp}
  }
  return(test = total < (val + gamma))
}

##########################################SVP###################################################
##########################################SVP###################################################
##########################################SVP###################################################
##########################################SVP###################################################

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
#' @param y A numeric vector representing the univariate signal to be segmented.
#' @param gamma A numeric value used as a threshold in the validation function and as a penalty for each segment.
#' @param test A function of the form `function(y, gamma)` returning TRUE if the segment is valid. Default is `valid_OP`.
#' @param pruning Logical indicating whether to apply pruning to the candidate set during the dynamic programming. Default is TRUE.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{changepoints}{Integer vector indicating the ending index of each segment (i.e., positions of changepoints).}
#'   \item{nb}{Integer vector of length \code{length(y)}. At each position \code{t}, it records the number of candidates tested.}
#'   \item{costQ}{Numeric vector of length \code{length(y)}. Quadratic cost value at each time step.}
#'   \item{R}{A matrix of dimension \code{(length(y)+1) x 3} containing, for each time step: 
#'     \code{Q} (cumulative cost), \code{K} (number of segments), and \code{s} (previous changepoint).}
#' }
#' 
#' @examples
#' y <- c(rnorm(50, 0, 1), rnorm(50, 5, 1))
#' result <- smallest_valid_partitioning_VR(y, gamma = 10, test = valid_OP)
#' print(result$changepoints)
#'
#' @export
smallest_valid_partitioning_VR <- function(y, gamma, test = valid_OP, pruning = TRUE)
{
  n <- length(y)
  # Initialisation de la matrice (Q, K, s)
  R <- matrix(Inf, nrow = n+1, ncol = 3)
  colnames(R) <- c("Q", "K", "s")
  
  # Condition initiale t=0
  R[1, ] <- c(-gamma, 0, 0)
  
  nb <- integer(n)    # nb de candidats examinés à chaque t
  costQ <- numeric(n) # Q_t final pour chaque t
  
  # Cumulative sum for optimized calculations
  cs_y  <- c(0, cumsum(y))
  cs_y2 <- c(0, cumsum(y^2))
  
  INDEX <- 0
  best_K_star <- 1
  
  for (t in 1:n)
  {
    best_Q <- Inf         # best cost for the current t
    best_K <- Inf         # best segment number
    best_s <- NA_integer_ # last changepoint

    for (s in INDEX)
    {
      #####
      ##### DANGER LINE. Remove if necessary
      #####
      if (!test(y[(s+1):t], gamma))
      {
        if (pruning) INDEX <- setdiff(INDEX, s) # we remove the invalid candidate
        next                                    # we stop here and we go to the next s
      }
      
      segment_cost <- (cs_y2[t+1]-cs_y2[s+1]) - (cs_y[t+1]-cs_y[s+1])^2 / (t-s)
      candidate_Q <- R[s+1, "Q"] + segment_cost
      candidate_K <- R[s+1, "K"] + 1
      
      if (candidate_K < best_K || (candidate_K == best_K && candidate_Q < best_Q))
      {
        best_Q <- candidate_Q
        costQ[t] <- best_Q
        best_K <- candidate_K
        best_s <- s
        best_K_star <- best_K
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
    s_prev <- R[t + 1, "s"]            # back to the previous s
    t <- s_prev                        # next loop
  }
  
  list(changepoints = unname(changepoints),
       nb           = nb,
       costQ        = costQ,
       R            = R)
}
