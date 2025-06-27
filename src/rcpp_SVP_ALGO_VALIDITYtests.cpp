#include <Rcpp.h>

using namespace Rcpp;

//' SVP Method using C++
//' 
//' @title SVP Method
//' @description This function implements the SVP (Smallest Valid Partitioning) method in Rcpp.
//' @param y A numeric vector representing the data to segment.
//' @param gamma A numeric value used as a threshold in the validation function and as a penalty for each segment.
//' @param test A function of the form `function(y, gamma)` returning TRUE if the segment is valid.
//' @param pruning Logical indicating whether to apply pruning to the candidate set during the dynamic programming. Default is TRUE.
//' @return A list with (1) the changepoint elements (each last index of each segment in \code{changepoints}), (2) a vector `\code{nb} saving the number of non-pruned elements at each iteration, (3) a vector \code{lastIndexSet} containing the non-pruned indices at the end of the algo and (4) a vector \code{costQ} saving the optimal cost at each time step.
//' @export
// [[Rcpp::export]]
List rcpp_smallest_valid_partitioning_VR(NumericVector y,
                                        double gamma,
                                        Function test,
                                        bool pruning = true)
{
  int n = y.size();
  
  // Cumulative sums
  NumericVector cs_y(n + 1);
  NumericVector cs_y2(n + 1);
  cs_y[0]  = 0;
  cs_y2[0] = 0;
  for (int i = 0; i < n; ++i)
  {
    cs_y[i + 1]  = cs_y[i] + y[i];
    cs_y2[i + 1] = cs_y2[i] + y[i] * y[i];
  }
  
  // R matrix: rows = t=0..n, cols = Q, K, s
  NumericMatrix Rmat(n + 1, 3);
  Rmat(0, 0) = 0;  // Q
  Rmat(0, 1) = 0;  // K
  Rmat(0, 2) = 0;  // s
  
  IntegerVector nb(n);
  NumericVector costQ(n);
  
  std::vector<int> INDEX;
  INDEX.reserve(n + 1);
  INDEX.push_back(0);
  
  // Main DP
  for (int t = 1; t <= n; ++t)
  {
    double best_Q = R_PosInf;
    int    best_K = INT_MAX;
    int    best_s = NA_INTEGER;
    
    // iterate over candidate changepoints
    for (size_t idx = 0; idx < INDEX.size(); ++idx)
    {
      int s = INDEX[idx];
      // segment [s, t-1]
      NumericVector segment = y[Range(s, t - 1)];
      bool ok = as<bool>(test(segment, gamma));
      if (!ok) {
        if (pruning) INDEX.erase(INDEX.begin() + idx);
        --idx;
        continue;
      }
      double seg_cost = (cs_y2[t] - cs_y2[s]) - std::pow(cs_y[t] - cs_y[s], 2) / (t - s);
      double cand_Q = Rmat(s, 0) + seg_cost;
      int    cand_K = Rmat(s, 1) + 1;
      
      if (cand_K < best_K || (cand_K == best_K && cand_Q < best_Q))
      {
        best_Q = cand_Q;
        best_K = cand_K;
        best_s = s;
      }
    }
    
    nb[t - 1]      = INDEX.size();
    costQ[t - 1]   = best_Q;
    Rmat(t, 0)     = best_Q;
    Rmat(t, 1)     = best_K;
    Rmat(t, 2)     = best_s;
    
    INDEX.push_back(t);
  }
  
  // Reconstruct changepoints
  std::vector<int> cp;
  cp.reserve(n);
  int t = n;
  while (t > 0)
  {
    cp.push_back(t);
    t = Rmat(t, 2);
  }
  std::reverse(cp.begin(), cp.end());
  
  return List::create(
    _["changepoints"] = IntegerVector(cp.begin(), cp.end()),
    _["nb"]           = nb,
    _["costQ"]        = costQ,
    _["R"]            = Rmat
  );
}
