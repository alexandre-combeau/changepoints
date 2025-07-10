##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' Optimal Partitioning algorithm using PELT
#' 
#' @title Optimal Partitioning using PELT
#' @description This function implements the OP algorithm using PELT of a given vector `data` with a given penalty term.
#' It finds the optimal changepoints that minimize the cost function using dynamic programming.
#' 
#' @param data A numeric vector representing the data to segment.
#' @param penalty A double value representing the penalty term for adding a new segment.
#' 
#' @return A list with (1) the changepoint elements (each last index of each segment in \code{changepoints}), (2) a vector `\code{nb} saving the number of non-pruned elements at each iteration, (3) a vector \code{lastIndexSet} containing the non-pruned indices at the end of the algo and (4) a vector \code{costQ} saving the optimal cost at each time step.
#' 
#' @examples
#' n <- 1000
#' data <- rep(c(0, 5, 2.5, 7), each = n) + rnorm(4 * n)
#' penalty <- 2 * log(length(data))
#' PELT_R <- pelt(data, penalty)
#' plot_segmentation(data, PELT_R$changepoints, title = "PELT Segmentation")
#' 
#' @export
pelt <- function(data, penalty)
{
  n <- length(data)
  
  # Initialize the costs and the changepoints
  Q <- rep(Inf, n + 1)
  Q[1] <- -penalty
  cp <- integer(n + 1)
  P <- 1
  length_P <- integer(n)
  
  # Cumulative sum for optimized calculations
  S1 <- c(0, cumsum(data))
  S2 <- c(0, cumsum(data^2))
  
  for (t in seq_len(n)) # avoid 1:0 issue
  {
    t1 <- t + 1
    costs <- rep(Inf, length(P))
    best_cost <- Inf
    
    for (i in seq_along(P))
    {
      s <- P[i]
      gaussian_cost <- (S2[t1] - S2[s]) - ((S1[t1] - S1[s])^2) / (t1 - s)
      costs[i] <- Q[s] + gaussian_cost + penalty
      if (costs[i] < best_cost)
      {
        best_cost <- costs[i]
        arg_min <- s
      }
    }
    
    Q[t1]       <- best_cost
    cp[t1]      <- arg_min
    length_P[t] <- length(P)
    
    # Pruning
    P <- c(P[costs <= Q[t1] + penalty], t1)
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
              lastIndexSet = rev(P)[-1],
              nb = length_P,
              costQ = Q[-1]))
}
