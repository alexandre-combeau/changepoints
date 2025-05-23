// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_optimal_partitioning
List rcpp_optimal_partitioning(NumericVector x, double beta);
RcppExport SEXP _changepoints_rcpp_optimal_partitioning(SEXP xSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_optimal_partitioning(x, beta));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_pelt
List rcpp_pelt(NumericVector x, double penalty, int minseglen);
RcppExport SEXP _changepoints_rcpp_pelt(SEXP xSEXP, SEXP penaltySEXP, SEXP minseglenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< int >::type minseglen(minseglenSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_pelt(x, penalty, minseglen));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_changepoints_rcpp_optimal_partitioning", (DL_FUNC) &_changepoints_rcpp_optimal_partitioning, 2},
    {"_changepoints_rcpp_pelt", (DL_FUNC) &_changepoints_rcpp_pelt, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_changepoints(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
