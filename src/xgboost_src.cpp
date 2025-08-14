#include <RcppArmadillo.h>
#include "shared_functions.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
List cox_custom_obj(const arma::vec& theta, 
                    const arma::vec& delta, 
                    const arma::vec& penalty_weights,
                    double lambda) {
  
  int n = delta.n_rows;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0 = rev_cumsum(exp_theta);
  arma::vec grad(n, arma::fill::zeros);
  arma::vec hess(n, arma::fill::zeros);
  
  // Gradient calculation
  arma::vec cum_risk_ratio = arma::zeros<arma::vec>(n);
  
  for(int i = 0; i < n; ++i){
    if(delta(i) == 1){
      cum_risk_ratio(i) = 1.0 / S0(i);
    }
  }
  
  arma::vec temp = rev_cumsum(cum_risk_ratio);
  grad = delta - exp_theta % temp;
  
  // Hessian (Diagonal only)
  arma::vec cum_risk_ratio_hess = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; ++i){
    if(delta(i) == 1){
      cum_risk_ratio_hess(i) = 1.0 / (S0(i) * S0(i));
    }
  }
  
  arma::vec temp_hess = rev_cumsum(cum_risk_ratio_hess);
  hess = - exp_theta % (temp - exp_theta % temp_hess);

// Optional penalty term added here (if penalty_weights provided):
  grad = grad - lambda * penalty_weights % theta;
  hess = hess - lambda * penalty_weights;
  
  return List::create(Named("grad") = grad,
                      Named("hess") = hess);
}