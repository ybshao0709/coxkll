#include <RcppArmadillo.h>
// #include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

arma::vec rev_cumsum(const arma::vec& X) {
  return arma::flipud(arma::cumsum(arma::flipud(X)));
}



/*
 * Computes the log-partial likelihood for the Cox-KL model.
 *
 * Parameters:
 *   lp              - Linear predictor (n × 1), i.e. Z * beta
 *   delta           - Event indicator vector (length n)
 *   delta_tilde     - Augmented failure indicator for KL divergence penalty
 *   eta             - Tuning parameter controlling KL divergence strength
 *   ind_start       - Starting indices of each stratum
 *   n_each_stratum  - Number of observations in each stratum
 *   beta            - Coefficient vector (p × 1)
 *   lambda          - Ridge Regularization parameter (≥ 0; default: 0.0)
 *
 * Returns:
 *   Scalar double: value of the stratified Cox-KL (potentially with ridge) log-partial likelihood.
 */
double coxkl_loglik(const arma::vec& lp,      
                    const arma::vec& delta,
                    const arma::vec& delta_tilde,
                    double eta,
                    const arma::uvec& ind_start,
                    const arma::vec& n_each_stratum,
                    const arma::vec& beta,
                    const double lambda = 0.0) {

  const arma::uword S = ind_start.n_elem;
  double loglik = 0.0;

  for (arma::uword s = 0; s < S; ++s) {
    const arma::uword start = ind_start(s);
    const arma::uword len   = static_cast<arma::uword>(n_each_stratum(s));
    if (len == 0) continue;
    const arma::uword end   = start + len - 1;

    arma::vec lp_s    = lp.subvec(start, end);
    arma::vec delta_s = delta.subvec(start, end);
    arma::vec dt_s    = delta_tilde.subvec(start, end);

    const double m = lp_s.max();
    arma::vec exp_lp_s_shift = arma::exp(lp_s - m);
    arma::vec S0_s = rev_cumsum(exp_lp_s_shift);

    for (arma::uword i = 0; i < len; ++i) {
      const double w_i = (delta_s(i) + eta * dt_s(i)) / (1.0 + eta);
      loglik += w_i * lp_s(i);
      if (delta_s(i) > 0) {
        loglik -= std::log(S0_s(i)) + m;
      }
    }
  }

  if (lambda > 0.0) {
    loglik -= 0.5 * lambda * arma::dot(beta, beta);
  }
  return loglik;
}


/*
 * Performs one Newton–Raphson update step for β in the Cox–KL model,
 * with optional backtracking line search.
 *
 *
 * Parameters:
 *   Z              - Covariate matrix (n × p)
 *   delta          - Event indicator vector (length n, 1 = event, 0 = censored)
 *   delta_tilde    - Augmented failure indicator used in KL weighting
 *   beta           - Current coefficient vector (p × 1)
 *   eta            - KL tuning parameter (η ≥ 0); η = 0 reduces to standard Cox
 *   ind_start      - Start indices of each stratum (length = #strata)
 *   n_each_stratum - Number of observations in each stratum
 *   lambda         - Ridge Regularization parameter (≥ 0; default: 0.0)
 *   backtrack      - If true, use a backtracking line search on the Newton step
 *
 * Returns:
 *   Newton step vector of step sizes.
 */
arma::vec BetaUpdate(const arma::mat& Z,
                     const arma::vec& delta,
                     const arma::vec& delta_tilde,
                     const arma::vec& beta,
                     const double eta,
                     const arma::uvec& ind_start,
                     const arma::vec& n_each_stratum,
                     const double lambda = 0.0,
                     bool backtrack = false) {

  const arma::uword p = Z.n_cols;       
  const arma::uword S = n_each_stratum.n_elem;

  arma::vec theta = Z * beta;
  // arma::vec exp_theta = arma::exp(theta);
  arma::vec L1(p, arma::fill::zeros);
  arma::mat L2(p, p, arma::fill::zeros);

  for (arma::uword s = 0; s < S; ++s) {
    const arma::uword start = ind_start(s);
    const arma::uword len = n_each_stratum(s);
    if (len == 0) continue;
    const arma::uword end = start + len - 1;

    arma::vec delta_s = delta.subvec(start, end);
    arma::vec delta_tilde_s = delta_tilde.subvec(start, end);
    arma::mat Z_s = Z.rows(start, end);

    arma::vec theta_s = theta.subvec(start, end);

    const double m = theta_s.max();
    arma::vec exp_theta_s = arma::exp(theta_s - m);
    arma::vec S0_s = rev_cumsum(exp_theta_s);

    arma::mat S1_s(len, p, arma::fill::zeros);
    arma::cube S2_s(len, p, p, arma::fill::zeros);

    for (arma::uword j = 0; j < p; ++j) {
      arma::vec Zj_s = Z_s.col(j);
      S1_s.col(j) = rev_cumsum(Zj_s % exp_theta_s);
      for (arma::uword k = 0; k < p; ++k) {
        arma::vec Zk_s = Z_s.col(k);
        S2_s.slice(j).col(k) = rev_cumsum(Zj_s % Zk_s % exp_theta_s);
      }
    }

    for (arma::uword i = 0; i < len; ++i) {
      const double d_i = delta_s(i);
      const double dt_i = delta_tilde_s(i);
      const double s0_i = S0_s(i);

      arma::rowvec z_i = Z_s.row(i);
      arma::rowvec s1_i = S1_s.row(i);

      L1 += ((d_i + eta * dt_i) / (1.0 + eta)) * z_i.t()
            - d_i * (s1_i / s0_i).t();

      if (d_i > 0) {
        arma::mat s2_i(p, p);
        for (arma::uword j = 0; j < p; ++j) {
          s2_i.col(j) = S2_s.slice(j).row(i).t();
        }
        arma::vec s1_div_s0 = s1_i.t() / s0_i;
        arma::mat outer = s1_div_s0 * s1_div_s0.t();
        L2 += d_i * (s2_i / s0_i - outer);
      }
    }
  }

  // Apply ridge to score & Hessian (no-op if lambda == 0)
  arma::vec L1_pen = L1 - lambda * beta;
  arma::mat L2_pen = L2;
  L2_pen.diag() += lambda;

  // Small numerical ridge (independent of lambda)
  L2_pen.diag() += 1e-8;

  arma::vec d_beta;
  bool ok = arma::solve(d_beta, L2_pen, L1_pen, arma::solve_opts::likely_sympd);

  if (!ok) {
    arma::mat L2_reg = L2_pen;
    L2_reg.diag() += 1e-4;
    ok = arma::solve(d_beta, L2_reg, L1_pen, arma::solve_opts::likely_sympd);
    if (!ok) Rcpp::stop("Hessian solve failed.");
  }

  if (backtrack) {
    const double armijo_s = 0.01;
    const double shrink_t = 0.60;
    const int max_bt = 10;
    
    double step_size = 1.0;

    double f0 = coxkl_loglik(theta, delta, delta_tilde, eta, ind_start, n_each_stratum, beta, lambda);
    double slope = arma::dot(L1_pen, d_beta);
    for (int bt = 0; bt < max_bt; ++bt) {
      arma::vec beta_try = beta + step_size * d_beta;
      arma::vec theta_try = Z * beta_try;

      double f_try = coxkl_loglik(theta_try, delta, delta_tilde, eta, ind_start, n_each_stratum, beta_try, lambda);
      if (f_try >= f0 +armijo_s * step_size * slope) break;
      step_size *= shrink_t;
    }
    d_beta *= step_size;
  }

  return d_beta;
}



 /*
 * Estimates Cox-KL model coefficients using Newton-Raphson iteration.
 *
 * Parameters:
 *   z              - Internal covariate matrix (n × p)
 *   delta          - Event indicator vector (length n), 1 = event, 0 = censored
 *   delta_tilde    - The predicted event indicator vector, defined using external risk scores (length n)
 *   n_each_stratum - Number of observations in each stratum (length S)
 *   eta            - Integration weight for KL penalty
 *   tol            - Convergence tolerance for Newton-Raphson updates (default: 1e-7)
 *   Mstop          - Maximum number of iterations (default: 50)
 *   lambda         - Regularization parameter (default: 0.0)
 *   backtrack      - Whether to use backtracking line search (default: false)
 * 
 * Returns:
 *   Estimated coefficient vector (β̂) for the Cox-KL model
 */
// [[Rcpp::export]]
arma::vec KL_Cox_Estimate_cpp(const arma::mat& Z, 
                              const arma::vec& delta, 
                              const arma::vec& delta_tilde,
                              const arma::vec& n_each_stratum,
                              const double eta,
                              arma::vec beta_initial, 
                              const double tol = 1.0e-7,
                              const int maxit = 50,
                              const double lambda = 0.0,
                              bool backtrack = false, 
                              bool message = false){
  arma::vec beta = beta_initial;
  
  const arma::uword S = n_each_stratum.n_elem;
  arma::uvec ind_start(S); 
  ind_start(0) = 0;
  for (arma::uword i = 1; i < S; ++i) {
    ind_start(i) = ind_start(i - 1) + n_each_stratum(i - 1);
  }

  for (int iter = 0; iter < maxit; ++iter) {
    arma::vec d_beta = BetaUpdate(Z, delta, delta_tilde, beta, eta, ind_start, n_each_stratum, lambda, backtrack);
    beta += d_beta;

    if (arma::max(arma::abs(d_beta)) < tol) {
      if (message){
        if (lambda == 0.0) {
          Rcpp::Rcout << "CoxKL converged after " << iter + 1 << " iterations." << std::endl;
        } else {
          Rcpp::Rcout << "CoxKL with ridge penalty (lambda=" << lambda << ") converged after " << iter + 1 << " iterations." << std::endl;
        }
      }
      break;
    }
  }
  return beta;
}


