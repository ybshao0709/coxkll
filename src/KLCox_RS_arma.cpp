#include <RcppArmadillo.h>
#include "shared_functions.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
double ddloglik_KL_RS_test(double z){
  return z;
}

// [[Rcpp::export]]
List loss_fn_cpp(const arma::mat& Z, const arma::vec& delta, arma::vec& beta) {
  int n = delta.n_rows;
  double loglik = 0;
  
  arma::vec theta = Z * beta;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0 = rev_cumsum(exp_theta);

  loglik = arma::accu(delta % (theta - arma::log(S0)));
  
  //print L2 dimensions:
  // Rcpp::Rcout << "L2 dimensions: " << L2.n_rows << " x " << L2.n_cols << std::endl;
  
  return List::create(Named("loglik") = loglik);
}


// [[Rcpp::export]]
List ddloglik_KL_RS_score(const arma::mat& Z, const arma::vec& delta, arma::vec& beta, const arma::vec& theta_tilde, const double &eta) {
  int p = beta.n_rows;
  int n = delta.n_rows;
  double loglik = 0;
  arma::vec L1 = arma::zeros<arma::vec>(p);
  arma::vec theta = Z * beta;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0 = rev_cumsum(exp_theta);
  arma::mat S1 = arma::zeros<arma::mat>(n, p);

  arma::vec exp_theta_tilde = arma::exp(theta_tilde);
  arma::vec S0_tilde = rev_cumsum(exp_theta_tilde);
  arma::mat S1_tilde = arma::zeros<arma::mat>(n, p);
  
  for (int i = 0; i < p; i++) {
    arma::vec S1_pre = Z.col(i) % exp_theta;
    S1.col(i) = rev_cumsum(S1_pre);
    
    arma::vec S1_pre_tilde = Z.col(i) % exp_theta_tilde;
    S1_tilde.col(i) = rev_cumsum(S1_pre_tilde);
  }
  
  // Perform the operation in steps for clarity
  arma::mat temp_1(n, p, arma::fill::zeros); 
  
  for (int i = 0; i < p; i++) {
    temp_1.col(i) = delta % (Z.col(i) + eta * (S1_tilde.col(i) / S0_tilde) - (1 + eta) * (S1.col(i) / S0));
  }
  
  L1 = arma::sum(temp_1, 0).t();
  loglik = arma::accu(delta % (theta - arma::log(S0)));
  
  //print L2 dimensions:
  // Rcpp::Rcout << "L2 dimensions: " << L2.n_rows << " x " << L2.n_cols << std::endl;
  
  return List::create(Named("loglik") = loglik,
                      Named("L1") = L1);
}

// // [[Rcpp::export]]
// List ddloglik_KL_RS_score_speed(const arma::mat& Z, const arma::vec& delta, arma::vec& beta, const arma::vec& theta_tilde, double eta) {
//   int p = beta.n_rows;
//   int n = delta.n_rows;
//   double loglik = 0;
//   
//   arma::vec theta = Z * beta;
//   
//   arma::vec exp_theta = arma::exp(theta);
//   arma::vec exp_theta_tilde = arma::exp(theta_tilde);
//   
//   arma::vec S0 = rev_cumsum(exp_theta);
//   arma::vec S0_tilde = rev_cumsum(exp_theta_tilde);
//   arma::mat S10       = Z % arma::repmat(exp_theta,1,p);
//   arma::mat S10_tilde = Z % arma::repmat(exp_theta_tilde,1,p);
//   
//   arma::vec Lambda    = arma::cumsum(delta/S0);
//   arma::mat matLambda = arma::repmat(Lambda,1,p);
//   arma::mat matdelta  = repmat(delta,1,p);
//   
//   arma::mat L         = matdelta%Z - S10%matLambda;
//   arma::vec L1        = arma::sum(L,0).t();
//   
//   loglik = arma::accu(delta%((theta)-log(S0)));
//   
// 
//   return List::create(Named("loglik") = loglik,
//                       Named("L1") = L1);
// }


// [[Rcpp::export]]
List ddloglik_KL_RS(const arma::mat& Z, const arma::vec& delta, arma::vec& beta, const arma::vec& theta_tilde, double &eta) {
  Rcpp::Rcout << "function started1" << std::endl;
  int p = beta.n_rows;
  int n = delta.n_rows;
  double loglik = 0;
  arma::vec L1 = arma::zeros<arma::vec>(p);
  arma::mat L2 = arma::zeros<arma::mat>(p, p);
  arma::vec theta = Z * beta;
  //print theta:
  // Rcpp::Rcout << "Theta: " << theta << std::endl;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0 = rev_cumsum(exp_theta);
  arma::mat S1 = arma::zeros<arma::mat>(n, p);
  arma::mat S2_i = arma::zeros<arma::mat>(n, p);

  arma::vec exp_theta_tilde = arma::exp(theta_tilde);
  arma::vec S0_tilde = rev_cumsum(exp_theta_tilde);
  arma::mat S1_tilde = arma::zeros<arma::mat>(n, p);
  
  if (Z.has_nan() || Z.has_inf()) {
    Rcpp::Rcout << "Error: Z contains NA or Inf values." << std::endl;
    return List();
  }
  
  // Check for NA or Inf values in beta
  if (beta.has_nan() || beta.has_inf()) {
    Rcpp::Rcout << "Error: beta contains NA or Inf values." << std::endl;
    return List();
  }
  
  if (exp_theta.has_nan() || exp_theta.has_inf()) {
    Rcpp::Rcout << "Error: exp_theta contains NA or Inf values." << std::endl;
    return List();
  }
  
  if (S0.has_nan() || S0.has_inf()) {
    Rcpp::Rcout << "Error: S0 contains NA or Inf values." << std::endl;
    return List();
  }
  
  if (exp_theta_tilde.has_nan() || exp_theta_tilde.has_inf()) {
    Rcpp::Rcout << "Error: exp_theta_tilde contains NA or Inf values." << std::endl;
    return List();
  }
  
  if (S0_tilde.has_nan() || S0_tilde.has_inf()) {
    Rcpp::Rcout << "Error: S0_tilde contains NA or Inf values." << std::endl;
    return List();
  }
  
  

  for (int i = 0; i < p; i++) {
    arma::vec S1_pre = Z.col(i) % exp_theta;
    S1.col(i) = rev_cumsum(S1_pre);
    
    if (S1.col(i).has_nan() || S1.col(i).has_inf()) {
      Rcpp::Rcout << "Error: S1 contains NA or Inf values in column " << i << std::endl;
      return List();
    }

    arma::vec S1_pre_tilde = Z.col(i) % exp_theta_tilde;
    S1_tilde.col(i) = rev_cumsum(S1_pre_tilde);
  }

  // First, ensure dimensions match
  if (delta.n_rows != Z.n_rows || S1_tilde.n_rows != Z.n_rows || S1.n_rows != Z.n_rows ||
      S0_tilde.n_rows != Z.n_rows || S0.n_rows != Z.n_rows) {
    stop("Dimension mismatch");
  }

  if (any(S0_tilde == 0) || any(S0 == 0)) {
    stop("Division by zero encountered in S0_tilde or S0");
  }

  // Perform the operation in steps for clarity
  arma::mat temp_1(n, p, arma::fill::zeros); 

  for (int i = 0; i < p; i++) {
    temp_1.col(i) = delta % (Z.col(i) + eta * (S1_tilde.col(i) / S0_tilde) - (1 + eta) * (S1.col(i) / S0));
  }

  L1 = arma::sum(temp_1, 0).t();
  loglik = arma::accu(delta % (theta - arma::log(S0)));
  ///

  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      arma::vec S2_i1 = (Z.col(j) % exp_theta) % Z.col(i);
      arma::vec S2_i_colj = rev_cumsum(S2_i1);
      S2_i.col(j) = S2_i_colj;
      
      if (S2_i.col(j).has_nan() || S2_i.col(j).has_inf()) {
        Rcpp::Rcout << "Error: S2_i contains NA or Inf values in column " << j << std::endl;
        return List();
      }

      arma::vec V1_i_colj = (S2_i.col(j) / S0) - (S1.col(j) % S1.col(i) / arma::square(S0));
      arma::vec V_i_colj = V1_i_colj % delta;
      L2(i, j) = arma::accu(V_i_colj) * (1 + eta);
    }
  }

  //print L2 dimensions:
  // Rcpp::Rcout << "L2 dimensions: " << L2.n_rows << " x " << L2.n_cols << std::endl;

  return List::create(Named("loglik") = loglik,
                      Named("L1") = L1,
                      Named("L2") = L2,
                      Named("S0") = S0);
}


arma::vec BetaUpdate(const arma::mat& Z, const arma::vec& delta, arma::vec& beta, const arma::vec& theta_tilde, double &eta) {
  int p = beta.n_rows;
  int n = delta.n_rows;
  double loglik = 0;
  arma::vec L1 = arma::zeros<arma::vec>(p);
  arma::mat L2 = arma::zeros<arma::mat>(p, p);
  arma::vec theta = Z * beta;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0 = rev_cumsum(exp_theta);
  arma::mat S1 = arma::zeros<arma::mat>(n, p);
  arma::mat S2_i = arma::zeros<arma::mat>(n, p);
  
  arma::vec exp_theta_tilde = arma::exp(theta_tilde);
  arma::vec S0_tilde = rev_cumsum(exp_theta_tilde);
  arma::mat S1_tilde = arma::zeros<arma::mat>(n, p);
  
  for (int i = 0; i < p; i++) {
    arma::vec S1_pre = Z.col(i) % exp_theta;
    S1.col(i) = rev_cumsum(S1_pre);
    
    arma::vec S1_pre_tilde = Z.col(i) % exp_theta_tilde;
    S1_tilde.col(i) = rev_cumsum(S1_pre_tilde);
  }
  
  // Perform the operation in steps for clarity
  arma::mat temp_1(n, p, arma::fill::zeros); 
  
  for (int i = 0; i < p; i++) {
    temp_1.col(i) = delta % (Z.col(i) + eta * (S1_tilde.col(i) / S0_tilde) - (1 + eta) * (S1.col(i) / S0));
  }
  
  L1 = arma::sum(temp_1, 0).t();
  loglik = arma::accu(delta % (theta - arma::log(S0)));
  ///
  
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      arma::vec S2_i1 = (Z.col(j) % exp_theta) % Z.col(i);
      arma::vec S2_i_colj = rev_cumsum(S2_i1);
      S2_i.col(j) = S2_i_colj;
      
      arma::vec V1_i_colj = (S2_i.col(j) / S0) - (S1.col(j) % S1.col(i) / arma::square(S0));
      arma::vec V_i_colj = V1_i_colj % delta;
      L2(i, j) = arma::accu(V_i_colj) * (1 + eta);
    }
  }
  
  //print L2 dimensions:
  // Rcpp::Rcout << "L2 dimensions: " << L2.n_rows << " x " << L2.n_cols << std::endl;
  arma::vec temp = solve(L2, L1, arma::solve_opts::fast);
  
  return temp;
}

// [[Rcpp::export]]
arma::vec KL_Cox_Estimate_cpp(const arma::mat& z, 
                          const arma::vec& delta, 
                          const arma::vec& time, 
                          const arma::vec& RS_internal, 
                          double eta, 
                          const double tol = 1.0e-7,
                          const bool returnBeta = false) {
  int p = z.n_cols;
  arma::vec beta = arma::zeros<arma::vec>(p);
  
  while (true) {
    arma::vec temp = BetaUpdate(z, delta, beta, RS_internal, eta);
    beta += temp;
    
    // Rcpp::Rcout << "Beta: " << beta << std::endl;

    if (max(abs(temp)) < tol) break;
  }
  
  return beta;
}

