#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' PELT Method using C++
//' 
//' @title PELT Method
//' @description This function implements the PELT (Pruned Exact Linear Time) method in Rcpp.
//' @param x A numeric vector representing the data to segment.
//' @param penalty A double value representing the penalty term for adding a new segment.
//' @param minseglen The minimum length of the segment.
//' @return A list with (1) the changepoint elements (each last index of each segment in \code{changepoints}), (2) a vector `\code{nb} saving the number of non-pruned elements at each iteration, (3) a vector \code{lastIndexSet} containing the non-pruned indices at the end of the algo and (4) a vector \code{costQ} saving the optimal cost at each time step.
//' @export
// [[Rcpp::export]]
List rcpp_pelt(NumericVector x, double penalty, int minseglen = 1)
{
  int n = x.size();
  NumericVector S1(n + 1, 0.0), S2(n + 1, 0.0);
  
  for (int i = 0; i < n; i++)
  {
    S1[i + 1] = S1[i] + x[i];
    S2[i + 1] = S2[i] + x[i] * x[i];
  }
  
  NumericVector Q(n + 1, R_PosInf);
  Q[0] = -penalty;
  IntegerVector last_cp(n + 1, 0);
  std::vector<int> R = {0};
  std::vector<int> length_R(n);
  
  for (int t = 0; t < n; t++)
  {
    int t1 = t + 1;
    double best_cost = R_PosInf;
    int best_index = -1;
    
    for (int s : R)
    {
      int len = t - s + 1;
      
      if (len >= minseglen)
      {
        double sum_x = S1[t1] - S1[s];
        double sum_x2 = S2[t1] - S2[s];
        double cost = sum_x2 - (sum_x * sum_x) / len;
        double total_cost = Q[s] + cost + penalty;
        if (total_cost < best_cost)
        {
          best_cost = total_cost;
          best_index = s;
        }
      }
    }
    
    Q[t1] = best_cost;
    last_cp[t1] = best_index;
    
    std::vector<int> newR;
    for (int s : R)
    {
      int len = t - s + 1;
      
      if (len >= minseglen)
      {
        double sum_x = S1[t1] - S1[s];
        double sum_x2 = S2[t1] - S2[s];
        double cost = sum_x2 - (sum_x * sum_x) / len;
        if (Q[s] + cost <= Q[t1] + penalty)
        {
          newR.push_back(s);
        }
      }
    }
    
    newR.push_back(t1);
    R = newR;
    length_R[t] = R.size();
  }
  
  // On retrouve les changepoints
  std::vector<int> changepoints;
  int i = n;
  while (last_cp[i] > 0)
  {
    changepoints.push_back(last_cp[i]);
    i = last_cp[i];
  }
  std::reverse(changepoints.begin(), changepoints.end());
  changepoints.push_back(n);
  std::reverse(R.begin(), R.end());
  
  return List::create(
    Named("changepoints") = wrap(changepoints),
    Named("lastIndexSet") = wrap(R),
    Named("nb") = wrap(length_R),
    Named("costQ") = wrap(Q));
}
