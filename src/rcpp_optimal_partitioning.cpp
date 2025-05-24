#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Optimal Partitioning Method using C++
//' 
//' @title Optimal Partitioning Method
//' @description This function computes the optimal partitioning of a given vector x with a given penalty term beta. It finds the optimal changepoints that minimize the cost function using dynamic programming.
//' @param x A numeric vector representing the data to segment.
//' @param beta A double value representing the penalty term for adding a new segment.
//' @return A list with (1) the changepoint elements (each last index of each segment in \code{changepoints}), (2) a vector `\code{nb} saving the number of non-pruned elements at each iteration, (3) a vector \code{lastIndexSet} containing the non-pruned indices at the end of the algo and (4) a vector \code{costQ} saving the optimal cost at each time step.
//' @export
// [[Rcpp::export]]
List rcpp_optimal_partitioning(NumericVector x, double beta)
{
  int n = x.size();
  
  // Initialisation
  std::vector<double> Q(n + 1, std::numeric_limits<double>::infinity());
  Q[0] = 0;
  std::vector<int> lastChange(n + 1, 0);
  std::vector<int> R(n);
  
  // Sommes cumulatives pour les calculs optimisés
  std::vector<double> cs_x(n + 1, 0.0);
  std::vector<double> cs_x2(n + 1, 0.0);
  
  for (int i = 0; i < n; i++)
  {
    R[i] = i + 1;
    cs_x[i + 1] = cs_x[i] + x[i];
    cs_x2[i + 1] = cs_x2[i] + x[i] * x[i];
  }
  
  // Calcul des coûts pour chaque sous-segment
  for (int t = 1; t <= n; t++)
  {
    for (int s = 0; s < t; s++)
    {
      // Coût de la moyenne quadratique sur le segment [s+1, t]
      double segment_cost = (cs_x2[t] - cs_x2[s]) -
        (cs_x[t] - cs_x[s]) * (cs_x[t] - cs_x[s]) / (t - s);
      
      // Coût total avec pénalité beta
      double cost = Q[s] + segment_cost + beta;
      
      // Mise à jour
      if (cost < Q[t])
      {
        Q[t] = cost;
        lastChange[t] = s;
      }
    }
  }
  
  // Reconstruction des points de rupture
  std::vector<int> optimal_cpts;
  int t = n;
  while (t > 0)
  {
    optimal_cpts.push_back(lastChange[t]);
    t = lastChange[t];
  }
  std::reverse(optimal_cpts.begin(), optimal_cpts.end());
  optimal_cpts.erase(optimal_cpts.begin());
  optimal_cpts.push_back(n);
  
  std::vector<int> last_R = std::vector<int>(R.rbegin(), R.rend());
  
  return List::create(
    Named("changepoints") = wrap(optimal_cpts),
    Named("lastIndexSet") = wrap(last_R),
    Named("nb") = wrap(R),
    Named("costQ") = wrap(Q));
}
