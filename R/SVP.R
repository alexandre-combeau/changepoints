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
smallest_valid_partitioning <- function(y, gamma)
{
  n <- length(y)
  
  # Initialisation de la matrice (Q,K,s)
  R <- matrix(Inf, nrow = n+1, ncol = 3)
  colnames(R) <- c("Q", "K", "s")
  
  # Condition initiale t=0
  R[1, ] <- c(0, 0, 0)
  
  nb <- integer(n) # nb de candidats examinés à chaque t
  costQ <- numeric(n) # Q_t final pour chaque t
  
  # Cumulative sum for optimized calculations
  cs_y  <- c(0, cumsum(y))
  cs_y2 <- c(0, cumsum(y^2))
  sums <- precompute_sums(y)
  
  for (t in 1:n)
  {
    best_Q <- Inf # meilleur coût trouvé pour ce t
    best_K <- Inf # meilleur nombre de segments
    best_s <- NA_integer_ # dernier changepoint retenu
    cnt <- 0 # nombre de candidats considérés
    
    for (s in 0:(t-1))
    {
      if (validity_function2_fast(y, s, t, gamma, sums = sums))
      {
        cnt <- cnt + 1
        segment_cost <- (cs_y2[t+1]-cs_y2[s+1]) - (cs_y[t+1]-cs_y[s+1])^2 / (t-s)
        
        candidate_Q <- R[s+1, "Q"] + segment_cost
        candidate_K <- R[s+1, "K"] + 1
        
        # on choisit d'abord le meilleur K, et en cas d'égalité le meilleur coût
        if (candidate_K < best_K || candidate_K == best_K && candidate_Q < best_Q)
        {
          best_Q <- candidate_Q
          best_K <- candidate_K
          best_s <- s
        }
      }
    }
    
    # on stocke le meilleur candidat pour R_t
    R[t+1, ] <- c(best_Q, best_K, best_s)
    
    nb[t] <- cnt
    costQ[t] <- best_Q
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
       costQ        = costQ,
       Matrix       = R)
}

# Example cost function : sum of squared deviations from mean
cost_function <- function(y, s, t)
{
  segment <- y[(s+1):t]
  mean_seg <- mean(segment)
  sum((segment - mean_seg)^2)
}


# Example validity function : standard deviation <= gamma
validity_function <- function(y, s, t, gamma)
{
  segment <- y[(s+1):t]
  if (length(segment) == 1)
  {
    sd_seg <- 0 # or an other value if we want to be strict
  }
  
  else
  {
    sd_seg <- sd(segment)
  }
  return(sd_seg <= gamma)
}

# Fonction pour pré-calculer les sommes cumulées
precompute_sums <- function(y) {
  list(cumsum_y  = c(0, cumsum(y)),
       cumsum_y2 = c(0, cumsum(y^2)))
}

# Fonction coût quadratique d’un segment [i:j]
fast_cost <- function(i, j, sums) {
  n <- j - i
  sum_y <- sums$cumsum_y[j + 1] - sums$cumsum_y[i + 1]
  sum_y2 <- sums$cumsum_y2[j + 1] - sums$cumsum_y2[i + 1]
  sum_y2 - (sum_y^2) / n
}

# Version rapide de validity_function2
validity_function2_fast <- function(y, s, t, gamma, sums = NULL)
{
  if ((t - s) < 2)
  {
    return(TRUE)
  }
  
  if (is.null(sums))
  {
    sums <- precompute_sums(y)
  }
  
  # Segment total
  total_cost <- fast_cost(s, t, sums)
  
  # Tous les u possibles
  candidates <- (s + 1):(t - 1)
  
  # Coûts gauche et droite pour tous les u
  costs_left <- fast_cost(s, candidates, sums)
  costs_right <- fast_cost(candidates, t, sums)
  
  # Gains
  gains <- total_cost - (costs_left + costs_right)
  
  return(max(gains) <= gamma)
}
