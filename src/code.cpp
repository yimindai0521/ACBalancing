#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//'@import RcppArmadillo
//'@export
// [[Rcpp::export]]
arma::colvec soft_c(const arma::colvec& beta, double lambda){
  // Your function code goes here
  double normvalue = arma::norm(beta, 2);
  if (normvalue > lambda){
    return((1 - lambda / normvalue) * beta);
  }else{
    return(0);
  }
}

//'@import RcppArmadillo
//'@export
// [[Rcpp::export]]
double objvalue(const arma::mat& X, const arma::colvec& beta, double lambda){
  // Your function code goes here
  double fobject = sum(exp(X * beta)) + lambda * arma::norm(beta, 2);
  return fobject;
}

//'@import RcppArmadillo
//'@export
// [[Rcpp::export]]
arma::colvec proximaloptim(const arma::mat& X, double rate, double lambda, arma::colvec& beta_start, double convergence, int iteration){
  // Your function code goes here
  int p = X.n_cols;
  int r = X.n_rows;
  int iterationtime = 0;
  double eps = 1000.0;
  arma::colvec beta;
  arma::colvec gradient;
  beta = beta_start;

  while (eps > convergence) {
    arma::colvec beta_old = beta;
    arma::colvec gradient = exp(X.t() * beta);
    arma::mat beta_old_mat(r, p);
    beta_old_mat.zeros();
    beta_old_mat.each_row() += gradient.t();
    arma::mat A(r, r, arma::fill::eye);
    arma::mat gradientmat = X;
    gradientmat.each_row() %= gradient.t();
    arma::mat hessian = arma::inv(X * gradientmat.t() + rate * A);
    arma::colvec vt = beta - hessian * X * gradient;
    beta = soft_c(vt, lambda);
    double fobj_old = objvalue(X.t(), beta_old, lambda);
    double fobj_new = objvalue(X.t(), beta, lambda);
    eps = std::abs(fobj_old - fobj_new);
    iterationtime++;
    if(iterationtime >= 1000){
      eps = 0;
    }
  }
  return beta;
}
