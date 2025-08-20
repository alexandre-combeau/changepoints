#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Smallest Valid Partitioning with Validation and Pruning using Rcpp
//' 
//' @title Smallest Valid Partitioning with Validation and Pruning
//' @description This function implements a dynamic programming approach to segment a univariate signal into the smallest number of valid segments, according to a user-defined validation function. Each segment must pass a validity test (e.g., based on variance, range, etc.). The algorithm minimizes a quadratic cost subject to this constraint.
//' 
//' @param data A numeric vector representing the univariate signal to be segmented.
//' @param gamma A numeric value used as a threshold in the validation function and as a penalty for each segment.
//' @param test A function of the form `function(data, gamma)` returning TRUE if the segment is valid. Default is `valid_OP`.
//' @param prune_if_unvalid Logical. If TRUE (default), the algorithm applies *segment-wise validation*:
//' at each time step, it tests whether the candidate segment \code{data[(s+1):t]} is valid using the
//' user-defined function \code{test}. If the segment fails the test, the candidate \code{s} is removed
//' (pruned) from the set of possible changepoints. This accelerates computation by avoiding invalid
//' segment extensions. If FALSE, the algorithm skips this validation and considers all candidate
//' segments without checking their validity (which can be faster but may return invalid segments).
//' 
//' @return A list with the following components :
//' \describe{
//'   \item{changepoints}{Integer vector indicating the ending index of each segment (i.e., positions of changepoints).}
//'   \item{nb}{Integer vector of length \code{length(data)}. At each position \code{t}, it records the number of candidates tested.}
//'   \item{costQ}{Numeric vector of length \code{length(data)}. Quadratic cost value at each time step.}
//'   \item{R}{A matrix of dimension \code{(length(data)+1) x 3} containing, for each time step :
//'     \describe{
//'       \item{Q}{cumulative cost}
//'       \item{K}{number of segments}
//'       \item{s}{previous changepoint}
//'     }
//'   }
//' }
//' 
//' @export
// [[Rcpp::export]]
List smallest_valid_partitioning_rcpp(std::vector<double> data,
                                      double gamma,
                                      Function test,
                                      bool prune_if_unvalid = true)
{
  size_t n = data.size();
  
  // Initialization of the R matrix
  NumericMatrix R(n + 1, 3); // Q, K, s
  std::fill(R.begin(), R.end(), std::numeric_limits<double>::infinity());
  R(0, 0) = 0.0; // Default value
  R(0, 1) = 0.0; // Default value
  R(0, 2) = 0.0; // Default value
  
  std::vector<size_t> nb(n);    // nb de candidats examinés à chaque t
  std::vector<double> costQ(n); // cost vector
  
  // Cumulative sum for optimized calculations
  std::vector<double> S1(n + 1, 0.0);
  std::vector<double> S2(n + 1, 0.0);
  for (size_t i = 0; i < n; ++i)
  {
    S1[i + 1] = S1[i] + data[i];
    S2[i + 1] = S2[i] + data[i] * data[i];
  }
  
  std::vector<size_t> INDEX = {0};
  
  // declare once outside the loop
  double best_Q;
  size_t best_K;
  size_t best_s;
  size_t s;
  bool valid;
  double candidate_Q;
  size_t candidate_K;
  
  for (size_t t = 1; t <= n; t++)
  {
    best_Q = std::numeric_limits<double>::infinity();
    best_K = std::numeric_limits<size_t>::max();
    best_s = 0;
    
    std::vector<size_t> new_INDEX;
    for (size_t k = 0; k < INDEX.size(); ++k)
    {
      s = INDEX[k];
      if (s >= t) continue;
      
      std::vector<double> seg(data.begin() + s, data.begin() + t);
      valid = as<bool>(test(seg, gamma));
      
      if (valid)
      {
        candidate_Q = R(s, 0) + (S2[t] - S2[s]) - (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s);
        candidate_K = R(s, 1) + 1;
        
        if (candidate_K < best_K || (candidate_K == best_K && candidate_Q < best_Q))
        {
          best_Q = candidate_Q;
          best_K = candidate_K;
          best_s = s;
          costQ[t - 1] = best_Q;
        }
        new_INDEX.push_back(s);
      }
    }
    
    // Pruning step
    std::vector<size_t> pruned_INDEX;
    
    for (size_t k = 0; k < new_INDEX.size(); ++k)
    {
      s = new_INDEX[k];
      if (s == t) continue;
      
      candidate_Q = R(s, 0) + (S2[t] - S2[s]) - (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s);
      candidate_K = R(s, 1);
      
      if (!(candidate_Q > best_Q && candidate_K == best_K))
      {
        pruned_INDEX.push_back(s);
      }
      
      else if (prune_if_unvalid)
      {
        pruned_INDEX.push_back(s);
      }
    }
    
    nb[t - 1] = INDEX.size();
    pruned_INDEX.push_back(t);
    INDEX = pruned_INDEX;
    
    R(t, 0) = best_Q;
    R(t, 1) = best_K;
    R(t, 2) = best_s;
  }
  
  // Reconstruct changepoints
  std::vector<size_t> changepoints;
  size_t i = n;
  while (i > 0)
  {
    changepoints.push_back(i);
    i = R(i, 2);
  }
  
  std::reverse(changepoints.begin(), changepoints.end());
  std::reverse(INDEX.begin(), INDEX.end());
  
  return List::create(
    _["changepoints"] = changepoints,
    _["lastIndexSet"] = INDEX,
    _["nb"] = nb,
    _["costQ"] = costQ,
    _["R"] = R(Range(1, R.nrow() - 1), Range(0, R.ncol() - 1))
  );
}
