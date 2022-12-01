// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// soft_c
arma::colvec soft_c(const arma::colvec& beta, double lambda);
RcppExport SEXP _ACBalancing_soft_c(SEXP betaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(soft_c(beta, lambda));
    return rcpp_result_gen;
END_RCPP
}
// objvalue
double objvalue(const arma::mat& X, const arma::colvec& beta, double lambda);
RcppExport SEXP _ACBalancing_objvalue(SEXP XSEXP, SEXP betaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(objvalue(X, beta, lambda));
    return rcpp_result_gen;
END_RCPP
}
// proximaloptim
arma::colvec proximaloptim(const arma::mat& X, double rate, double lambda, arma::colvec& beta_start, double convergence, int iteration);
RcppExport SEXP _ACBalancing_proximaloptim(SEXP XSEXP, SEXP rateSEXP, SEXP lambdaSEXP, SEXP beta_startSEXP, SEXP convergenceSEXP, SEXP iterationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta_start(beta_startSEXP);
    Rcpp::traits::input_parameter< double >::type convergence(convergenceSEXP);
    Rcpp::traits::input_parameter< int >::type iteration(iterationSEXP);
    rcpp_result_gen = Rcpp::wrap(proximaloptim(X, rate, lambda, beta_start, convergence, iteration));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ACBalancing_soft_c", (DL_FUNC) &_ACBalancing_soft_c, 2},
    {"_ACBalancing_objvalue", (DL_FUNC) &_ACBalancing_objvalue, 3},
    {"_ACBalancing_proximaloptim", (DL_FUNC) &_ACBalancing_proximaloptim, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_ACBalancing(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}