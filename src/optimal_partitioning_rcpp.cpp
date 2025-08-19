#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Optimal Partitioning algorithm using Rcpp
//' 
//' @title Optimal Partitioning
//' @description This function implements the OP algorithm of a given vector `data` with a given penalty term.
//' It finds the optimal changepoints that minimize the cost function using dynamic programming.
//' 
//' @param data A numeric vector representing the data to segment.
//' @param penalty A double value representing the penalty term for adding a new segment.
//' 
//' @return A list with (1) the `\code{changepoints}` elements, (2) a vector `\code{nb}` saving the number of non-pruned elements at each iteration, (3) a vector `\code{lastIndexSet}` containing the non-pruned indices at the end of the algo and (4) a vector `\code{costQ}` saving the optimal cost at each time step.
//' 
//' @examples
//' n <- 1000
//' data <- rep(c(0, 5, 2.5, 7), each = n) + rnorm(4 * n)
//' penalty <- 2 * log(length(data))
//' OP_Rcpp <- optimal_partitioning_rcpp(data, penalty)
//' plot_segmentation(data, OP_Rcpp$changepoints, title = "OP Segmentation")
//' 
//' @export
// [[Rcpp::export]]
List optimal_partitioning_rcpp(std::vector<double> data, double penalty)
{
  size_t n = data.size();
  
  // Initialize the costs and the changepoints
  std::vector<double> Q(n + 1, std::numeric_limits<double>::infinity());
  Q[0] = -penalty;
  std::vector<size_t> lastChange(n + 1, 0);
  std::vector<size_t> R(n);
  
  // Cumulative sum for optimized calculations
  std::vector<double> S1(n + 1, 0.0), S2(n + 1, 0.0);
  for (size_t i = 0; i < n; i++)
  {
    R[i] = i + 1;
    S1[i + 1] = S1[i] + data[i];
    S2[i + 1] = S2[i] + data[i] * data[i];
  }
  
  // Cost calculation for each sub-segment
  
  double tempQ;
  
  for (size_t t = 1; t <= n; t++)
  {
    for (size_t s = 0; s < t; s++)
    {
      // Total cost with beta penalisation
      tempQ = Q[s] + (S2[t] - S2[s]) - (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s) + penalty;
      
      // Minimization
      if (tempQ < Q[t])
      {
        Q[t] = tempQ;
        lastChange[t] = s;
      }
    }
  }
  
  // Changepoints reconstruction (backtracking)
  std::vector<size_t> optimal_cpts;
  size_t t = n;
  while (t > 0)
  {
    optimal_cpts.push_back(lastChange[t]);
    t = lastChange[t];
  }
  std::reverse(optimal_cpts.begin(), optimal_cpts.end());
  optimal_cpts.erase(optimal_cpts.begin());
  optimal_cpts.push_back(n);
  
  std::vector<size_t> last_R(R.rbegin(), R.rend());
  
  return List::create(
    Named("changepoints") = optimal_cpts,
    Named("lastIndexSet") = last_R,
    Named("nb")           = R,
    Named("costQ")        = std::vector<double>(Q.begin() + 1, Q.end()));
}
