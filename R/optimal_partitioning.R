##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' Optimal Partitioning Method using R
#'
#' @title Optimal Partitioning Method
#' @description This function computes the optimal partitioning of a given vector x with a given penalty term beta.
#' It finds the optimal changepoints that minimize the cost function using dynamic programming.
#' 
#' @param data A numeric vector representing the data to segment.
#' @param beta A double value representing the penalty term for adding a new segment.
#' 
#' @return A list with (1) the changepoint elements (each last index of each segment in \code{changepoints}), (2) a vector `\code{nb} saving the number of non-pruned elements at each iteration, (3) a vector \code{lastIndexSet} containing the non-pruned indices at the end of the algo and (4) a vector \code{costQ} saving the optimal cost at each time step.
#'
#' @examples
#' n <- 1000
#' data <- rep(c(0, 5, 2.5, 7), each = n) + rnorm(4 * n)
#' beta <- 2 * log(length(data))
#' OP_R <- optimal_partitioning(data, beta)
#' plot_segmentation(data, OP_R$changepoints, title = "OP Segmentation")
#' 
#' @export
optimal_partitioning <- function(data, beta)
{
  n <- length(data)
  
  # Initialize the costs and the changepoints
  Q <- rep(Inf, n + 1)
  Q[1] <- -beta
  lastChange <- rep(0, n + 1)
  
  # Cumulative sum for optimized calculations
  cs_x  <- c(0, cumsum(data))
  cs_x2 <- c(0, cumsum(data^2))
  
  # Cost calculation for each sub-segment
  for (t in 1:n)
  {
    for (s in 0:(t - 1))
    {
      # Segment cost [s+1, t]
      segment_cost <- (cs_x2[t+1]-cs_x2[s+1]) - (cs_x[t+1]-cs_x[s+1])^2 / (t-s)
      
      # Total cost with beta penalisation
      cost <- Q[s + 1] + segment_cost + beta
      
      # Minimization
      if (cost < Q[t + 1])
      {
        Q[t + 1] <- cost
        lastChange[t + 1] <- s
      }
    }
  }
  
  # Changepoints reconstruction (backtracking)
  optimal_cpts <- integer(0)
  t <- n
  while (t > 0)
  {
    optimal_cpts <- c(lastChange[t + 1], optimal_cpts)
    t <- lastChange[t + 1]
  }
  
  changepoints = c(optimal_cpts[-1], n)
  
  return(list(changepoints = changepoints,
              lastIndexSet = n:1,
              nb = 1:n,
              costQ = Q[-1]))
}
