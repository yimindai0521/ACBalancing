#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec soft_c(const arma::colvec& beta, double lambda){ // proximal operator for l_2 norm
  // Your function code goes here
  double normvalue = arma::norm(beta, 2);
  if (normvalue > lambda){
    return((1 - lambda / normvalue) * beta); // output = (1 - lambda / ||beta||_2) * beta if ||beta||_2 >= lambda
  }else{
    return(0); // output = 0 if lambda > ||beta||_2
  }
}

// [[Rcpp::export]]
double objvalue(const arma::mat& X, const arma::colvec& beta, double lambda){ // Calculate the value of loss function
  // Your function code goes here
  double fobject = sum(exp(X * beta)) + lambda * arma::norm(beta, 2); // loss = exp(X * beta) + lambda * ||beta||_2
  return fobject;
}

// [[Rcpp::export]]
arma::colvec proximaloptim(const arma::mat& X, double rate, double lambda, arma::colvec& beta_start, double convergence, int iteration){
  // Your function code goes here
  int p = X.n_cols; // the number of columns in X
  int r = X.n_rows; // the number of rows in X
  int iterationtime = 0; // initialize maximum iteration time
  double eps = 1000.0; // initialize the absolute tolerance used to determine a stopping criteri
  arma::colvec beta; // initialize beta
  arma::colvec gradient; // initialize gradient
  beta = beta_start; // initialize beta_start

  while (eps > convergence) {
    arma::colvec beta_old = beta; // Save beta in beta_start
    arma::colvec gradient = exp(X.t() * beta); // Calculate gradient
    arma::mat beta_old_mat(r, p); // initialize beta_old matrix
    beta_old_mat.zeros(); // initialize beta_old matrix to be zero
    beta_old_mat.each_row() += gradient.t(); // Set each row of beta_old matrix to be gradient
    arma::mat A(r, r, arma::fill::eye); //  initialize identity matrix
    arma::mat gradientmat = X; // initialize gradient matrix to be X
    gradientmat.each_row() %= gradient.t(); // each row of gradient matrix to time gradient
    arma::mat hessian = arma::inv(X * gradientmat.t() + rate * A); // Tikhonov-regularized Newton update
    arma::colvec vt = beta - hessian * X * gradient; // update beta using gradient descent
    beta = soft_c(vt, lambda); // proximal operator to find new beta
    double fobj_old = objvalue(X.t(), beta_old, lambda); // Calculate the value of loss function with beta_old
    double fobj_new = objvalue(X.t(), beta, lambda); // Calculate the value of loss function with beta
    eps = std::abs(fobj_old - fobj_new); // Calculate the difference of loss function after update beta
    iterationtime++; // iteration time = iteration time + 1
    if(iterationtime >= iteration){
      eps = 0; // if iterationtime is larger than maximum time, stop iteration.
    }
  }
  return beta;
}
