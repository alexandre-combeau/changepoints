#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Optimal Partitioning algorithm using PELT
//' 
//' @title Optimal Partitioning using PELT
//' @description This function implements the OP algorithm using PELT of a given vector `data` with a given penalty term.
//' It finds the optimal changepoints that minimize the cost function using dynamic programming.
//' 
//' @param data A numeric vector representing the data to segment.
//' @param penalty A double value representing the penalty term for adding a new segment.
//' 
//' @return A list with
//' (1) the changepoint elements (each last index of each segment in \code{changepoints}),
//' (2) a vector \code{nb} saving the number of non-pruned elements at each iteration,
//' (3) a vector \code{lastIndexSet} containing the non-pruned indices at the end of the algorithm,
//' (4) a vector \code{costQ} saving the optimal cost at each time step.
//' 
//' @examples
//' n <- 1000
//' data <- rep(c(0, 5, 2.5, 7), each = n) + rnorm(4 * n)
//' penalty <- 2 * log(length(data))
//' PELT_Rcpp <- pelt_rcpp(data, penalty)
//' plot_segmentation(data, PELT_Rcpp$changepoints, title = "PELT Segmentation")
//' 
//' @export
// [[Rcpp::export]]
List pelt_rcpp(std::vector<double> data, double penalty)
{
  int n = data.size();
  
  // Initialize the costs and the changepoints
  std::vector<double> Q(n + 1, std::numeric_limits<double>::infinity());
  Q[0] = -penalty;
  std::vector<int> last_cp(n + 1, 0);
  std::vector<int> P(1, 0);
  std::vector<int> length_P(n);
  
  // Cumulative sum for optimized calculations
  std::vector<double> S1(n + 1, 0.0), S2(n + 1, 0.0);
  for (int i = 0; i < n; i++)
  {
    S1[i + 1] = S1[i] + data[i];
    S2[i + 1] = S2[i] + data[i] * data[i];
  }
  
  for (int t = 0; t < n; t++)
  {
    int t1 = t + 1;
    double best_cost = R_PosInf;
    std::vector<double> costs(P.size(), std::numeric_limits<double>::infinity());
    int arg_min = -1;
    
    for (size_t i = 0; i < P.size(); i++)
    {
      int s = P[i];
      double sum_x = S1[t1] - S1[s];
      double sum_x2 = S2[t1] - S2[s];
      double gaussian_cost = sum_x2 - (sum_x * sum_x) / (t1 - s);
      costs[i] = Q[s] + gaussian_cost + penalty;
      if (costs[i] < best_cost)
      {
        best_cost = costs[i];
        arg_min = s;
      }
    }
    
    Q[t1] = best_cost;
    last_cp[t1] = arg_min;
    length_P[t] = P.size();
    
    // Pruning
    std::vector<int> newP;
    for (size_t i = 0; i < P.size(); i++)
    {
      int s = P[i];
      if (costs[i] <= Q[t1] + penalty) newP.push_back(s);
    }
    newP.push_back(t1);
    P = newP;
  }
  
  // Backtracking
  std::vector<int> changepoints;
  int i = n;
  while (last_cp[i] > 0)
  {
    changepoints.push_back(last_cp[i]);
    i = last_cp[i];
  }
  std::reverse(changepoints.begin(), changepoints.end());
  changepoints.push_back(n);
  std::reverse(P.begin(), P.end());
  
  return List::create(
    Named("changepoints") = wrap(changepoints),
    Named("lastIndexSet") = wrap(P),
    Named("nb") = wrap(length_P),
    Named("costQ") = wrap(std::vector<double>(Q.begin() + 1, Q.end())));
}
