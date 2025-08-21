#include "tests.h"
#include <Rcpp.h>
#include <functional>
using namespace Rcpp;

// [[Rcpp::export]]
List smallest_valid_partitioning_rcpp_seq_old(NumericVector data,
                                              double gamma,
                                              std::string test,
                                              bool all_full_validity = true)
{
  int n = data.size();

  NumericMatrix R(n + 1, 3); // Q, K, s
  std::fill(R.begin(), R.end(), R_PosInf);
  R(0, 0) = 0.0;
  R(0, 1) = 0.0;
  R(0, 2) = 0.0;

  IntegerVector nb(n);
  NumericVector costQ(n);

  // Precompute cumulative sums
  NumericVector cs_y(n + 1, 0.0);
  NumericVector cs_y2(n + 1, 0.0);
  for (int i = 0; i < n; ++i) {
    cs_y[i + 1] = cs_y[i] + data[i];
    cs_y2[i + 1] = cs_y2[i] + data[i] * data[i];
  }

  // Precompute segment costs: seg_cost[t][s] for 0 <= s < t <= n
  std::vector<std::vector<double>> seg_cost(n + 1, std::vector<double>(n + 1, 0.0));
  for (int t = 1; t <= n; ++t) {
    for (int s = 0; s < t; ++s) {
      double sum_y = cs_y[t] - cs_y[s];
      double sum_y2 = cs_y2[t] - cs_y2[s];
      seg_cost[t][s] = sum_y2 - (sum_y * sum_y) / (t - s);
    }
  }

  // Define a generic initializer for test objects
  std::function<std::unique_ptr<TestBase>()> newTest;
  if (test == "gaussian_mean") {
    newTest = []() { return std::make_unique<GaussianMean>(); };
  }
  // Add more cases here for other tests, e.g.:
  // else if (test == "bernoulli_mean") {
  //   newTest = []() { return std::make_unique<BernoulliMean>(); };
  // }
  else {
    stop("Unknown test type");
  }

  std::vector<int> INDEX = {0};
  std::vector<std::unique_ptr<TestBase>> TESTS;
  TESTS.push_back(newTest());

  for (int t = 1; t <= n; t++) {
    double best_Q = R_PosInf;
    int best_K = INT_MAX;
    int best_s = -1;

    std::vector<int> new_INDEX;
    std::vector<std::unique_ptr<TestBase>> new_TESTS;

    for (size_t k = 0; k < INDEX.size(); ++k) {
      int s = INDEX[k];
      auto& test_instance = TESTS[k];

      if (s >= t) continue;

      test_instance->update(data[t - 1]);
      double stat_val = test_instance->statistic();
      bool valid = true;
      if (all_full_validity) {
        valid = stat_val < gamma;
      }

      if (valid) {
        double candidate_Q = R(s, 0) + seg_cost[t][s];
        int candidate_K = R(s, 1) + 1;

        if (candidate_K < best_K || (candidate_K == best_K && candidate_Q < best_Q)) {
          best_Q = candidate_Q;
          best_K = candidate_K;
          best_s = s;
          costQ[t - 1] = best_Q;
        }
        new_INDEX.push_back(s);
        new_TESTS.push_back(std::move(test_instance));
      }
    }

    // Pruning step
    std::vector<int> pruned_INDEX;
    std::vector<std::unique_ptr<TestBase>> pruned_TESTS;
    for (size_t k = 0; k < new_INDEX.size(); ++k) {
      int s = new_INDEX[k];
      auto& test_instance = new_TESTS[k];
      if (s == t) continue;

      double candidate_Q = R(s, 0) + seg_cost[t][s];
      int candidate_K = R(s, 1);

      if (!(candidate_Q > best_Q && candidate_K == best_K)) {
        pruned_INDEX.push_back(s);
        pruned_TESTS.push_back(std::move(test_instance));
      }
    }
    // Add new candidate for t
    pruned_INDEX.push_back(t);
    pruned_TESTS.push_back(newTest());
    
    nb[t - 1] = INDEX.size();
    INDEX = pruned_INDEX;
    TESTS = std::move(pruned_TESTS);

    R(t, 0) = best_Q;
    R(t, 1) = best_K;
    R(t, 2) = best_s;
  }

  // Reconstruct changepoints
  IntegerVector changepoints;
  int t = n;
  while (t > 0) {
    changepoints.push_front(t);
    t = R(t, 2);
  }

  return List::create(
    _["changepoints"] = changepoints,
    _["nb"] = nb,
    _["costQ"] = costQ,
    _["R"] = R(Range(1, R.nrow() - 1), Range(0, R.ncol() - 1))
  );
}
