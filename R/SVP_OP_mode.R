##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' Smallest Valid Partitioning Method using R
#' 
#' @title SVP Method
#' @description This function implements the SVP method in R.
#' It computes the optimal partitioning of a given numeric vector `x` using a penalized cost function.
#' @param y A numeric vector representing the data to segment.
#' @param gamma Threshold: segment is valid if validity_function(.) lower than gamma
#' @return list(
#'    changepoints = last index of each segment (in order),
#'    nb           = number of non-pruned candidates at every t = 1..n,
#'    costQ        = optimal cumulative cost Q_t, t = 1..n
#' )
#' 
#' @export
smallest_valid_partitioning_op_mode <- function(y, gamma)
{
  n <- length(y)
  
  # Initialisation de la matrice (Q,K,s)
  R <- matrix(Inf, nrow = n+1, ncol = 3)
  colnames(R) <- c("Q", "K", "s")
  
  # Condition initiale t=0
  R[1, ] <- c(-gamma, 0, 0)
  
  nb <- integer(n) # nb de candidats examinés à chaque t
  costQ <- numeric(n) # Q_t final pour chaque t
  
  # Cumulative sum for optimized calculations
  cs_y  <- c(0, cumsum(y))
  cs_y2 <- c(0, cumsum(y^2))
  
  for (t in 1:n)
  {
    best_Q <- Inf # meilleur coût trouvé pour ce t
    best_K <- Inf # meilleur nombre de segments
    best_s <- NA_integer_ # dernier changepoint retenu
    cnt <- 0 # nombre de candidats considérés
    
    for (s in 0:(t-1))
    {
      cnt <- cnt + 1
      segment_cost <- (cs_y2[t+1]-cs_y2[s+1]) - (cs_y[t+1]-cs_y[s+1])^2 / (t-s)
      
      candidate_Q <- R[s+1, "Q"] + segment_cost + gamma
      candidate_K <- R[s+1, "K"] + 1
      
      # on choisit d'abord le meilleur K, et en cas d'égalité le meilleur coût
      if (candidate_Q < best_Q)
      {
        best_Q <- candidate_Q
        costQ[t] <- best_Q
        best_K <- candidate_K
        best_s <- s
      }
    }
    
    # on stocke le meilleur candidat pour R_t
    R[t+1, ] <- c(best_Q, best_K, best_s)
    
    nb[t] <- cnt
  }
  
  ## Reconstruction
  changepoints <- integer(0)
  t <- n
  while (t > 0)
  {
    changepoints <- c(t, changepoints) # dernier indice du segment courant
    s_prev <- R[t + 1, "s"] # retour au précédent s
    t <- s_prev # prochaine boucle
  }
  
  list(changepoints = changepoints,
       nb           = nb,
       costQ        = costQ)
}