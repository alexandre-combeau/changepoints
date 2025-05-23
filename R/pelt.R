##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' PELT Method using R
#' 
#' @title PELT Method
#' @description This function implements the PELT (Pruned Exact Linear Time) method in pure R.
#' It computes the optimal partitioning of a given numeric vector `x` using a penalized cost function.
#' @param x A numeric vector representing the data to segment.
#' @param penalty A double value representing the penalty term for adding a new segment.
#' @param minseglen The minimum length of the segment
#' @return An integer vector of indices representing the optimal breakpoints.
#' @examples
#' x <- c(1, 2, 3, 10, 10, 10, 20, 20, 20, 5, 5, 5)
#' penalty <- 5
#' pelt(x, penalty)
#' 
#' @export
pelt <- function(x, penalty, minseglen = 1) {
  n <- length(x)
  
  # Cumulative sum for optimized calculations
  S1 <- c(0, cumsum(x))
  S2 <- c(0, cumsum(x^2))
  
  Q <- rep(Inf, n + 1)
  total_cost <- Inf
  Q[1] <- -penalty
  cp <- integer(n + 1)
  R <- 1
  length_R <- integer(n)
  
  for (t in seq_len(n)) # seq_len pour eviter les problemes comme 1:0
  {
    t1 <- t + 1
    costs <- rep(Inf, length(R))
    
    for (i in seq_along(R))
    {
      s <- R[i]
      len <- t - s + 1
      if (len >= minseglen)
      {
        gaussian_cost <- (S2[t1] - S2[s]) - ((S1[t1] - S1[s])^2) / len
        costs[i] <- Q[s] + gaussian_cost + penalty
        cost_best <- costs[i]
        if (cost_best < total_cost)
        {
          total_cost <- cost_best
          arg_min <- s
        }
      }
    }
    
    Q[t1]       <- costs[arg_min]
    cp[t1]      <- R[arg_min]
    length_R[t] <- length(R)
    
    # Pruning
    R <- c(R[costs <= Q[t1] + penalty], t1)
  }
  
  # Backtracking
  changepoints <- integer(0)
  i <- n + 1
  while (cp[i] > 1)
  {
    changepoints <- c(cp[i] - 1, changepoints)
    i <- cp[i]
  }
  
  changepoints = c(changepoints, n)
  return(list(changepoints = changepoints,
              lastIndexSet = rev(R),
              nb = length_R,
              costQ = Q))
}
