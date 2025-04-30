// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
#include <RcppArmadilloExtensions/sample.h>
#include "shared_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;


double norm(const arma::vec &x, const int &p) {
  double x_norm = 0;
  for (int j=0; j<p; j++) {
    x_norm = x_norm + pow(x(j),2);
  }
  x_norm = sqrt(x_norm);
  return(x_norm);
}

// Cross product of the jth column of x with y
double crossprod(const arma::mat &x, const arma::vec &y, const int &n, const int &j) {
  double val = 0;
  for (int i=0; i<n; i++) val += x(i,j) * y(i);
  return(val);
}

// [[Rcpp::export]]
double  maxgrad(const arma::mat &X, const arma::vec &y, const arma::vec &K, const arma::vec &m){
  
  //intialize
  int n = X.n_rows;
  int J = K.n_elem - 1;
  double zmax = 0;

  for (int g = 0; g < J; ++g)
  {
    int Kg = K(g+1) - K(g);
    arma::vec Z = arma::zeros<arma::vec>(Kg);
    for (int j = K(g); j < K(g+1); ++j)
    { 
      Z(j - K(g)) = crossprod(X, y, n, j);
    }
    double z = norm(Z, Kg) / m(g);
    if ( z > zmax) zmax = z;
  }
  
  return zmax;
}

// Soft-thresholding operator
double S(double z, double l) {
  if (z > l) return(z-l);
  if (z < -l) return(z+l);
  return(0);
}


// standardize function
// [[Rcpp::export]]
List standardize(const arma::mat &x) {
	int n = x.n_rows;
	int p = x.n_cols;
	arma::mat XX = arma::zeros<arma::mat>(n,p) ;
	arma::vec c = arma::zeros<arma::vec>(p), s = arma::zeros<arma::vec>(p);

	for (int j = 0; j < p; ++j)
	{
		//center:
		c(j) = 0;
		for (int i = 0; i < n; ++i)
		{
			c(j) = c(j) + x(i,j);
		}
		c(j) = c(j)/n;
		for (int i = 0; i < n; ++i)
		{
			XX(i,j) = x(i,j) - c(j);
		}

		//Scale:
		s(j) = 0;
		for (int i = 0; i < n; ++i)
		{
			s(j) += pow(XX(i,j),2);
		}
		s(j) = sqrt(s(j)/n);
		for (int i = 0; i < n; ++i)
		{
			XX(i,j) = XX(i,j)/s(j);
		}
	}

  return List::create(_["XX"]=XX,
                      _["c"]=c,
                      _["s"]=s);

}




// Group descent update -- cox
// call in main loop: gd_cox(beta, X, r, eta, v, g, K1, n, l, p, penalty, l1, l2, gamma, df, a, maxChange, actNum);
// r: first order derivative (n*1)
// eta: X*beta
// v: =1, hessian matrix approximation
// g: index of groups
// l: lambda
//	 lam1 = lambda(l) * group_multiplier(g) * alpha;
//	 lam2 = 0, (lambda(l) * group_multiplier(g) * (1-alpha);)
void gd_cox(arma::mat &b, const arma::mat &x, arma::vec &r, arma::vec &eta, double &v, const int &g,
            const arma::uvec &K1, const int &n, int &l, int &p, const std::string &penalty,  double &lam1,
            double &lam2, const double &gamma, arma::vec &df, arma::vec &a, double &maxChange) {

  // Calculate z
  int K = K1(g+1) - K1(g);        //In group g, there are K covaraites
  arma::vec z = arma::zeros<arma::vec>(K);          //z_g = 1/n * (X_g^T * r) + beta_g; a K-dimensional vector.
  for (int j=K1(g); j<K1(g+1); j++) {z(j-K1(g)) = crossprod(x, r, n, j)/n + a(j);}  
  double z_norm = norm(z,K);      //||z_g||_2,  L2 norm.
  
  // Calculate X_g^T * S (where S is the score vector r) for the group
  arma::vec XgTr = arma::zeros<arma::vec>(K);
  for (int j = K1(g); j < K1(g+1); j++) {
    // Use arma::dot, assuming 'r' is the score vector S passed from gdfit_cox_kl
    XgTr(j - K1(g)) = arma::dot(x.col(j), r);
  }
  
  // Update b based on the penalty
  // if (penalty == "lasso") {
  // Lasso penalty (L1)
  double len = S(v * z_norm, lam1) / (v * (1 + lam2));
  if (len != 0 || a(K1(g)) != 0) {
    for (int j = K1(g); j < K1(g+1); j++) {
      b(l, j) = len * z(j - K1(g)) / z_norm;
      double shift = b(l, j) - a(j);
      if (fabs(shift) > maxChange) maxChange = fabs(shift);
      for (int i = 0; i < n; i++) {
        double si = shift * x(i, j);
        r(i) -= si;
        eta(i) += si;
      }
    }
  }
  if (len > 0) df(l) += K * len / z_norm;
  // } else if (penalty == "ridge") {
  //   // Ridge penalty (L2)
  //   // Rcout << "Ridge penalty" << std::endl;
  //   // Rcout << "lam2: " << lam2 << std::endl;
  //   for (int j = K1(g); j < K1(g+1); j++) {
  //     b(l, j) = (a(j) + z(j - K1(g))) / (v * (1 + lam2));
  //     double shift = b(l, j) - a(j);
  //     if (fabs(shift) > maxChange) maxChange = fabs(shift);
  //     for (int i = 0; i < n; i++) {
  //       double si = shift * x(i, j);
  //       r(i) -= si;
  //       eta(i) += si/n;
  //     }
  //   }
  //   df(l) += K;  // Degrees of freedom for ridge
  //   // double denominator_ridge = v + lam2;
  //   // 
  //   // if (denominator_ridge <= 1e-15) {
  //   //   Rcpp::warning("Ridge update denominator near zero for group %d. Check v (%f) and lam2 (%f). Setting update to zero.", g, v, lam2);
  //   //   // Set beta for this group to zero or keep old value? Zeroing might be safer.
  //   //   b.row(l).subvec(K1(g), K1(g+1)-1).zeros();
  //   // } else {
  //   //   arma::vec beta_new_g = (v * a.subvec(K1(g), K1(g+1)-1) + XgTr) / denominator_ridge;
  //   //   
  //   //   // Calculate shift and update eta, r, maxChange
  //   //   for (int j_idx = 0; j_idx < K; ++j_idx) {
  //   //     int j = K1(g) + j_idx; // Actual column index
  //   //     double shift = beta_new_g(j_idx) - a(j);
  //   //     b(l, j) = beta_new_g(j_idx); // Update beta matrix
  //   //     
  //   //     if (!std::isfinite(shift)) {
  //   //       Rcpp::warning("Non-finite shift in Ridge update for group %d, feature %d.", g, j);
  //   //       shift = 0;
  //   //     }
  //   //     if (fabs(shift) > maxChange) maxChange = fabs(shift);
  //   //     
  //   //     // Update eta and approximate update to r (score)
  //   //     for (int i = 0; i < n; i++) {
  //   //       double si = shift * x(i, j);
  //   //       eta(i) += si;
  //   //       r(i) -= si * n; // Approximate score update (using v)
  //   //     }
  //   //   }
  //   // }
  //   // df(l) += K; // Degrees of freedom for ridge
  // } else {
  //   // Handle other penalties if necessary
  //   Rcpp::stop("Unknown penalty type.");
  // }
  
  // Update b
  // double len;
  // len = S(v * z_norm, lam1) / (v * (1 + lam2));   //||z_g||_2 - lambda(l), rightnow v = 1, lam2 =0. so the denominator has no effect
  
  // if (len != 0 || a(K1(g)) != 0) {
  //   // If necessary, update the K-dim beta, r, eta;
  //   for (int j=K1(g); j<K1(g+1); j++) { 
  //     b(l,j) = len * z(j-K1(g)) / z_norm;
  //     double shift = b(l,j)-a(j);
  //     if (fabs(shift) > maxChange) maxChange = fabs(shift);
  //     for (int i=0; i<n; i++) {
  //       double si = shift*x(i,j);
  //       r(i)   	 -= si;
  //       eta(i) 	 += si;
  //     }
  //   }
  // }

  // Update df
  // if (len > 0) df(l) += K * len / z_norm;
}



// // [[Rcpp::export]]
// List gdfit_cox(const arma::mat &X, const arma::vec &d, const std::string &penalty,
//                 const arma::uvec &K1, const int &K0, const arma::vec &lambda, const double &alpha, const double &eps,
//                 const int &max_iter, const double &gamma, const arma::vec &group_multiplier, const int &dfmax, const int &gmax,
//                 const bool &warn, const bool &user,
//                 const int &actIter = 50){
// 
//   // cout<<"done 1"<<endl;
//   // Lengths/dimensions
//   int n = d.n_elem;                     //number of samples
//   int L = lambda.n_elem;                //number of penalized coefficients (for grid search)
//   int J = K1.n_elem - 1;                //number of groups
//   int p = X.n_cols;                     //number of predictors
// 
//   int tot_iter = 0;                     // iteration index
// 
//   // outcomes:
//   List res;
//   arma::mat beta = arma::zeros<arma::mat>(L,p);           // store L*p beta
//   arma::vec Loss = arma::zeros<arma::vec>(L);             // loss function for each lambda
//   arma::uvec iter = arma::zeros<arma::uvec>(L);             // iterations taken for each lambda
//   arma::vec df   = arma::zeros<arma::vec>(L);
//   arma::mat Eta  = arma::zeros<arma::mat>(L,n);
// 
//   // interrmediate quantities  add notes:
//   arma::vec a    = arma::zeros<arma::vec>(p);
//   arma::vec r    = arma::zeros<arma::vec>(n), h = arma::zeros<arma::vec>(n), haz = arma::zeros<arma::vec>(n), rsk = arma::zeros<arma::vec>(n), eta = arma::zeros<arma::vec>(n);
//   arma::vec e    = arma::zeros<arma::vec>(J);
// 
//   int lstart, ng, nv, violations;
//   double shift, l1, l2, nullDev, v, s, maxChange;
// 
//   // Initialization
//   // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
//   rsk(n-1) = 1;                            //rsk: (e.g.  n = 5,  rsk = (5,4,3,2,1)  )
//   for (int i = n-2; i >= 0; i--) {
//     rsk(i) = rsk(i+1) + 1;
//   }
//   nullDev = 0;
//   for (int i = 0; i < n; i++){
//     nullDev -= d(i)*log(rsk(i));           //nullDev : - sum(delta_i * log(rsk(i)))
//   }
//   if (user) {
//     lstart = 0;
//   } else {
//     lstart = 1;
//     Loss(0) = nullDev;                                     //Loss initialize as nullDev
//   }
// 
//   // cout<<"done 2"<<endl;
//   // Path through different  lambda:
//   for (int l = lstart; l < L; ++l){
// 
//     tot_iter = 0; // reset max_iter for each lambda
// 
//     Rcpp::checkUserInterrupt();
//     if(l != 0) {
//       // Assign previous beta to a (b: store all the beta)
//       for (int j = 0; j < p; ++j) a(j) = beta(l-1,j);
//     }
// 
//     // Check dfmax, gmax
//     ng = 0;
//     nv = 0;             // sotre all the number of covariates in each group
//     for (int g = 0; g < J; g++) {
//       if (a(K1(g)) != 0) {
//         ng++;
//         nv = nv + (K1(g+1)-K1(g));
//       }
//     }
//     if (ng > gmax || nv > dfmax) {
//       for (int ll=l; ll<L; ll++) iter(ll) = NA_INTEGER;
//       break;
//     }
// 
//     //begin iteration:
//     while(tot_iter < max_iter){
// 
//       //start to update the inner set
//       for (int hot_i = 0; hot_i < actIter; ++hot_i){
// 	      iter(l)++;
// 	      tot_iter++;
// 	      Loss(l) = 0;
// 	      df(l)   = 0;
// 
// 	      // Calculate haz, risk
// 	      haz = exp(eta);                                        // eta (to be updated below), for (int i=0; i<n; i++) haz(i) = exp(eta(i));
// 	      rsk(n-1) = haz(n-1);
// 	      for (int i=n-2; i>=0; i--) {
// 	        rsk(i) = rsk(i+1) + haz(i);                          // rsk : reverse culmulative sum of haz
// 	      }
// 	      for (int i=0; i<n; i++) {
// 	        Loss(l) += d(i)*eta(i) - d(i)*log(rsk(i));           // update Loss
// 	      }
// 
// 	      //Approximate L:
// 	      h(0) = d(0)/rsk(0);                                    // delta_0 / rsk0
// 	      v = 1;                                                 // identitical matrix to approximate
// 	      for (int i=1; i<n; i++) {
// 	        h(i) = h(i-1) + d(i)/rsk(i);                         // hi: sum(delta_l / rsk_l) (l: at risk set)   this is cumulative baseline hazard
// 	      }
// 	      for (int i=0; i<n; i++) {
// 	        h(i) = h(i)*haz(i);
// 	        s    = d(i) - h(i);                  // corresponds to our first order derivative wrt pseudo-response
// 	        if (h(i)==0) r(i)=0;                 // corresponds to our l1/eta derivative
// 	        else r(i) = s/v;                     // use identical matrix to approximate the hessian matrix
// 	      }
// 
// 	      // Check for saturation    //not sure if we should change this
// 	      if (Loss(l)/nullDev < .01) {
// 	        if (warn) warning("Model saturated; exiting...");
// 	        for (int ll=l; ll<L; ll++) iter(ll) = NA_INTEGER;
// 	        tot_iter = max_iter;
// 	        break;
// 	      }
// 
// 	      // Update unpenalized covariates
// 	      maxChange = 0;
// 	      for (int j=0; j<K0; j++) {
// 	        shift = crossprod(X, r, n, j)/n;
// 	        if (fabs(shift) > maxChange) maxChange = fabs(shift);
// 	        beta(l,j) = shift + a(j);
// 	        for (int i=0; i<n; i++) {
// 	          double si = shift * X(i,j);
// 	          r(i)   -= si;              // minus si ()
// 	          eta(i) += si;
// 	        }
// 	        df(l)++;
// 	      }
// 	      // Update penalized groups
// 	      for (int g=0; g<J; g++) {
// 	        if (penalty == "ridge") {
// 	          l1 = 0;
// 	          l2 = lambda(l) * group_multiplier(g);  // Ridge penalty uses lam2
// 	        } else if (penalty == "lasso") {
// 	          l1 = lambda(l) * group_multiplier(g);  // Lasso penalty uses lam1
// 	          l2 = 0;
// 	        } else {
// 	          l1 = lambda(l) * group_multiplier(g) * alpha;
// 	          l2 = lambda(l) * group_multiplier(g) * (1 - alpha);
// 	        }
// 	        if (e(g)!=0) gd_cox(beta, X, r, eta, v, g, K1, n, l, p, penalty, l1, l2,
// 	                         gamma, df, a, maxChange);
// 	      }
// 	      // Check convergence
// 	      for (int j=0; j<p; j++) a(j) = beta(l,j);
// 	      if (maxChange < eps) break;
//       }
// 
// 
//       //outer loop:
// 
//       // Scan for violations
//       violations = 0;
//       for (int g=0; g<J; g++) {
//         if (e(g)==0) {
//           // l1 = lambda(l) * group_multiplier(g) * alpha;
//           // l2 = lambda(l) * group_multiplier(g) * (1-alpha);
//           if (penalty == "ridge") {
//             l1 = 0;
//             l2 = lambda(l) * group_multiplier(g);  // Ridge penalty uses lam2
//           } else if (penalty == "lasso") {
//             l1 = lambda(l) * group_multiplier(g);  // Lasso penalty uses lam1
//             l2 = 0;
//           } else {
//             l1 = lambda(l) * group_multiplier(g) * alpha;
//             l2 = lambda(l) * group_multiplier(g) * (1 - alpha);
//           }
//           gd_cox(beta, X, r, eta, v, g, K1, n, l, p, penalty, l1, l2, gamma, df, a, maxChange);
//           if (beta(l,K1(g)) != 0) {
//             e(g) = 1;
//             violations++;
//           }
//         }
//       }
//       if (violations==0) {
//         for (int i=0; i<n; i++) Eta(l,i) = eta(i);  // can use .row
//         break;
//       }
//       for (int j=0; j<p; j++) a(j) = beta(l,j);     // can use .row
//     }
// 
//     if(tot_iter == max_iter){
//       Rcout <<"lambda "<< lambda(l) << " failed to converge within "<<max_iter<<" iterations"<<endl;
//     }
// 
//   }
// 
//   res["beta"]  = beta;
//   res["iter"]  = iter;
//   res["df"]    = df;
//   res["Loss"]  = Loss;
//   res["Eta"]   = Eta;
// 
//   return res;
// 
// }


// [[Rcpp::export]]
List gdfit_cox_kl(const arma::mat &X, const arma::vec &d,
                  const arma::vec &delta_tilde,
                  const arma::uvec &K1, const int &K0, const arma::vec &lambda, const double &alpha, const double &eps,
                  const double &eta_kl,
                  const int &max_iter, const double &gamma, const arma::vec &group_multiplier, const int &dfmax, const int &gmax,
                  const bool &warn, const bool &user,
                  const int &actIter = 20){

  // Rcpp::Rcout<<"done 1"<<endl;
  // Lengths/dimensions
  int n = d.n_elem;                     //number of samples
  int L = lambda.n_elem;                //number of penalized coefficients (for grid search)
  int J = K1.n_elem - 1;                //number of groups
  int p = X.n_cols;                     //number of predictors

  int tot_iter = 0;                     // iteration index

  // outcomes:
  List res;
  arma::mat beta = arma::zeros<arma::mat>(L,p);           // store L*p beta
  arma::vec Loss = arma::zeros<arma::vec>(L);             // loss function for each lambda
  arma::uvec iter = arma::zeros<arma::uvec>(L);             // iterations taken for each lambda
  arma::vec df   = arma::zeros<arma::vec>(L);
  arma::mat Eta  = arma::zeros<arma::mat>(L,n);

  // interrmediate quantities  add notes:
  arma::vec a    = arma::zeros<arma::vec>(p);
  arma::vec r    = arma::zeros<arma::vec>(n), h = arma::zeros<arma::vec>(n), haz = arma::zeros<arma::vec>(n), rsk = arma::zeros<arma::vec>(n), eta = arma::zeros<arma::vec>(n);
  arma::vec e    = arma::zeros<arma::vec>(J);

  int lstart, ng, nv, violations;
  double shift, l1, l2, nullDev, v, s, maxChange;

  // Initialization
  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  rsk(n-1) = 1;                            //rsk: (e.g.  n = 5,  rsk = (5,4,3,2,1)  )
  for (int i = n-2; i >= 0; i--) {
    rsk(i) = rsk(i+1) + 1;
  }
  nullDev = 0;
  for (int i = 0; i < n; i++){
    nullDev -= d(i)*log(rsk(i));           //nullDev : - sum(delta_i * log(rsk(i)))
  }
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
    Loss(0) = nullDev;                                     //Loss initialize as nullDev
  }

  // cout<<"done 2"<<endl;
  // Path through different  lambda:
  // Rcpp::Rcout<<"done 2"<<endl;
  // Rcpp::Rcout<<"lstart: "<<lstart<<endl;
  // Rcpp::Rcout<<"L: "<<L<<endl;
  for (int l = lstart; l < L; ++l){

    tot_iter = 0; // reset max_iter for each lambda

    Rcpp::checkUserInterrupt();
    if(l != 0) {
      // Assign previous beta to a (b: store all the beta)
      for (int j = 0; j < p; ++j) a(j) = beta(l-1,j);
    }

    // Check dfmax, gmax
    ng = 0;
    nv = 0;             // sotre all the number of covariates in each group
    for (int g = 0; g < J; g++) {
      if (a(K1(g)) != 0) {
        ng++;
        nv = nv + (K1(g+1)-K1(g));
      }
    }
    if (ng > gmax || nv > dfmax) {
      for (int ll=l; ll<L; ll++) iter(ll) = NA_INTEGER;
      break;
    }
    
    std::string penalty_str;
    if (alpha == 0.0) {
      Rcpp::warning("alpha should be larger than 0; defaulting to 'lasso'");
    }
    
    // Rcpp::Rcout << "begin: " << std::endl;
    //begin iteration:
    while(tot_iter < max_iter){

      //start to update the inner set
      for (int hot_i = 0; hot_i < actIter; ++hot_i){
	      iter(l)++;
	      tot_iter++;
	      Loss(l) = 0;
	      df(l)   = 0;


        //////////////////////////////////////////

	      // Calculate haz, risk
        // Rcpp::Rcout << "Loss: " << Loss(l) << std::endl;
        // Rcpp::Rcout << "r: " << r(0) << std::endl;
        // calculateRiskAndUpdateLoss(eta, d, Loss(l), r);
        // Rcpp::Rcout << "Loss: " << Loss(l) << std::endl;
        // Rcpp::Rcout << "r: " << r(0) << std::endl;
        // List res_haz_risk = calculateRiskAndUpdateLoss(eta, d, Loss(l), r);
        // Loss(l) = res_haz_risk["Loss"];
        // r = Rcpp::as<arma::vec>(res_haz_risk["r"]);
        // arma::vec r = res_haz_risk["r"];
        //
        // Calculate haz, risk


	      haz = exp(eta);                                        // eta (to be updated below), for (int i=0; i<n; i++) haz(i) = exp(eta(i));
	      // Rcpp::Rcout << "beta: " << beta << std::endl;
	      
	      rsk(n-1) = haz(n-1);
	      for (int i=n-2; i>=0; i--) {
	        rsk(i) = rsk(i+1) + haz(i);                          // Cumulative sum from the end: sum_{k=i}^{n-1} exp(z_k*beta) = sum_{k=i}^{n-1} haz_k
	      }
	      for (int i=0; i<n; i++) {
	        double weighted_term = (d(i)+eta_kl*delta_tilde(i))/(1.0+eta_kl); // weighted term
	        Loss(l) += weighted_term*eta(i) - d(i)*log(rsk(i));           // update Loss
	      }

	      //Approximate L:
	      h(0) = d(0)/rsk(0);                                    // delta_0 / rsk0
	      v = 1;                                                 // identitical matrix to approximate
	      for (int i=1; i<n; i++) {
	        h(i) = h(i-1) + d(i)/rsk(i);                         // hi: sum(delta_l / rsk_l) (l: at risk set)   this is cumulative baseline hazard
	      }
	      // Rcpp::Rcout << "h(i): " << h << std::endl;
	      for (int i=0; i<n; i++) {
	        h(i) = h(i)*haz(i);
	        s    = (d(i)+eta_kl*delta_tilde(i))/(1+eta_kl) - h(i);   // corresponds to our first order derivative wrt pseudo-response
	        if (h(i)==0) {
            // raiseerror if d(i) is not 0:
            if(d(i) != 0){
              // Rcpp::Rcout << "h(i): " << h(i) << std::endl;
              // Rcpp::Rcout << "d(i): " << d(i) << std::endl;
              // Rcpp::Rcout << "i: " << i << std::endl;
              Rcpp::stop("h(i)==0 but d(i) is not 0");
            }
            r(i) = (d(i)+eta_kl*delta_tilde(i))/(1+eta_kl);
          }                                     // corresponds to our l1/eta derivative
	        else r(i) = s/v;                      // use identical matrix to approximate the hessian matrix
	      }
        // Rcpp::Rcout << "Loss: " << Loss(l) << std::endl;
        // Rcpp::Rcout << "r: " << accu(r) << std::endl;
        //////////////////////////////////////////

	      // Check for saturation    //not sure if we should change this
	      if (!std::isfinite(Loss(l))) {
	        if (warn) warning("Loss is non-finite at lambda %f; exiting...", lambda(l));
	        for (int ll=l; ll<L; ll++) iter(ll) = NA_INTEGER;
	        tot_iter = max_iter; // Force exit from outer while loop
	        break; // Exit hot_i loop
	      }
	      if (Loss(l)/nullDev < .01) {
	        if (warn) warning("Model saturated; exiting...");
	        for (int ll=l; ll<L; ll++) iter(ll) = NA_INTEGER;
	        tot_iter = max_iter;
	        break;
	      }

	      // Update unpenalized covariates
	      maxChange = 0;
	      for (int j=0; j<K0; j++) {
	        shift = crossprod(X, r, n, j)/n;
	        if (fabs(shift) > maxChange) maxChange = fabs(shift);
	        beta(l,j) = shift + a(j);
	        for (int i=0; i<n; i++) {
	          double si = shift * X(i,j);
	          r(i)   -= si;              // minus si ()
	          eta(i) += si;
	        }
	        df(l)++;
	      }
	      // Update penalized groups
	      for (int g=0; g<J; g++) {
	        l1 = lambda(l) * group_multiplier(g) * alpha;
	        l2 = lambda(l) * group_multiplier(g) * (1-alpha);
	        if (e(g)!=0) gd_cox(beta, X, r, eta, v, g, K1, n, l, p, penalty_str, l1, l2,
	                         gamma, df, a, maxChange);
	      }
	      // Check convergence
	      for (int j=0; j<p; j++) a(j) = beta(l,j);
	      if (maxChange < eps) break;
      }


      //outer loop:

      // Scan for violations
      violations = 0;
      for (int g=0; g<J; g++) {
        if (e(g)==0) {
          l1 = lambda(l) * group_multiplier(g) * alpha;
          l2 = lambda(l) * group_multiplier(g) * (1-alpha);
          gd_cox(beta, X, r, eta, v, g, K1, n, l, p, penalty_str, l1, l2, gamma, df, a, maxChange);
          if (beta(l,K1(g)) != 0) {
            e(g) = 1;
            violations++;
          }
        }
      }
      if (violations==0) {
        for (int i=0; i<n; i++) Eta(l,i) = eta(i);  // can use .row
        break;
      }
      for (int j=0; j<p; j++) a(j) = beta(l,j);     // can use .row
    }

    if(tot_iter == max_iter){
      Rcout <<"lambda "<< lambda(l) << " failed to converge within "<<max_iter<<" iterations"<<endl;
    }
    
    // Rcpp::Rcout << "tot_iter: " << tot_iter << std::endl;
  }

  res["beta"]  = beta;
  res["iter"]  = iter;
  res["df"]    = df;
  res["Loss"]  = Loss;
  res["Eta"]   = Eta;

  return res;

}