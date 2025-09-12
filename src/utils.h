// shared_functions.h
#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

arma::vec rev_cumsum(const arma::vec& X);

arma::vec calculateDeltaTilde(const arma::vec& event, 
                              const arma::vec& time, 
                              const arma::vec& theta_tilde,
                              const arma::vec& n_each_stratum);

double pl_cal_theta(const arma::vec& lp,
                    const arma::vec& delta,
                    const arma::vec& time,
                    const arma::vec& n_each_stratum);

#endif // UTILS_H
