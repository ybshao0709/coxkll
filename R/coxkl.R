KL_Cox_Estimate <- function(z, delta, time, RS_internal, eta, tol=1.0e-7){
  
  p <- ncol(z)
  beta = as.matrix(rep(0,p))
  RS_internal = as.matrix(RS_internal)
  
  repeat{
    diff=ddloglik_KL_RS(z,delta, beta, RS_internal, eta)
    G=diff$L1
    H=diff$L2
    S0=diff$S0
    
    Lambda=cumsum(delta/S0)
    S=exp(-Lambda)
    
    temp=solve(H)%*%G
    beta=beta+temp
    
    if(max(abs(temp))<tol) break
  }
  
  return(beta)
}


#' Cox Proportional Hazards Model with KL Divergence for Data Integration
#'
#' Estimates the coefficients of a Cox proportional hazards model incorporating external information using Kullback-Leibler (KL) divergence. External information can be integrated via precomputed external risk scores (RS) or externally derived coefficient estimates (beta). The strength of integration is controlled by tuning parameters (eta).
#'
#' @param z Numeric covariate matrix with rows representing observations and columns representing predictor variables. All variables must be numeric.
#' @param delta Numeric event indicator vector (1 = event occurred, 0 = censored).
#' @param time Numeric vector of observed event or censoring times. No sorting is required.
#' @param RS Numeric vector or matrix of precomputed external risk scores. Length (or number of rows) must match the number of observations. If not provided, `beta` must be specified.
#' @param beta Numeric vector of externally derived coefficients (e.g., from prior studies). Length must match the number of columns in `z`. Use zeros to represent absent covariates if fewer coefficients are externally available. If not provided, `RS` must be specified.
#' @param eta_list Numeric vector of tuning parameters controlling the integration strength of the external information. Higher values indicate greater reliance on external information.
#' @param tol Numeric scalar controlling convergence tolerance for optimization. Default is `1.0e-7`.
#' @param Mstop Integer specifying the maximum number of iterations allowed for the optimization algorithm. Default is `50`.
#'
#' 
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{beta_list}}{List of estimated coefficient vectors corresponding to each `eta` value.}
#'     \item{\code{LP_list}}{List of linear predictor vectors computed on the training data corresponding to each `eta` value.}
#'     \item{\code{eta_list}}{Original vector of tuning parameters (`eta`) provided by the user.}
#'   }
#'
#' @examples
#' # Example data has 6 covariates: z1, z2, z5, z6 (continuous) and z3, z4 (binary)
#' # External beta is simulated from a homogeneous data distribution.
#'
#' # Using external beta information
#' result1 <- coxkl(z = ExampleData$z,
#'                  delta = ExampleData$status,
#'                  time = ExampleData$time,
#'                  beta = ExampleData$beta_external,
#'                  eta_list = seq(0, 5, 1))
#'
#' # Using external risk score information
#' rs <- as.matrix(ExampleData$z) %*% as.matrix(ExampleData$beta_external)
#' result2 <- coxkl(z = ExampleData$z,
#'                  delta = ExampleData$status,
#'                  time = ExampleData$time,
#'                  RS = rs,
#'                  eta_list = seq(0, 5, 1))
#'
#' @export

coxkl <- function(z, delta, time, RS = NULL, beta = NULL, eta_list, tol=1.0e-7, Mstop = 50){
  
  if(is.null(RS) && is.null(beta)) {
    stop("Error: No external information is provided. Either RS or beta must be provided.")
  }  else if(is.null(RS) && !is.null(beta)) {
    # Check if the dimension of beta matches the number of columns in z
    if(length(beta) == ncol(z)) {
      print("External beta information is used.")
      RS <- as.matrix(z) %*% as.matrix(beta)
    } else {
      stop("Error: The dimension of beta does not match the number of columns in z.")
    }
  } else if(!is.null(RS)) {
    print("External Risk Score information is used.")
  }
  
  time_order <- order(time)
  delta      <- delta[time_order]
  z          <- z[time_order,]
  RS <- RS[time_order]
  time       <- time[time_order]
  
  z_mat <- as.matrix(z)
  delta_mat <- as.matrix(delta)
  p     <- ncol(z_mat)
  beta  <- as.matrix(rep(0,p))
  RS    <- as.matrix(RS)
  
  beta_list <- list()
  LP_list <- list()
  for (i in seq_along(eta_list)){
    eta <- eta_list[i]
    beta_train <- KL_Cox_Estimate(z_mat, delta_mat, t_train, RS, eta=eta)
    LP_train <- as.matrix(z) %*% as.matrix(beta_train)
    
    beta_list[[i]] <- beta_train
    LP_list[[i]] <- LP_train
  }
  
  results <- list(beta_list = beta_list, LP_list = LP_list, eta_list = eta_list)
  return(results)
}



# coxkl with sorted data as input
coxkl_sorted <- function(z, delta, time, RS, eta_list, tol=1.0e-7, Mstop = 50){
  z_mat <- as.matrix(z)
  delta_mat <- as.matrix(delta)
  p     <- ncol(z_mat)
  beta  <- as.matrix(rep(0,p))
  RS    <- as.matrix(RS)
  
  beta_list <- list()
  LP_list <- list()
  
  for (i in seq_along(eta_list)){
    eta <- eta_list[i]
    beta_train <- KL_Cox_Estimate(z_mat, delta_mat, t_train, RS, eta=eta)
    LP_train <- as.matrix(z) %*% as.matrix(beta_train)
    
    beta_list[[i]] <- beta_train
    LP_list[[i]] <- LP_train
  }
  
  results <- list(beta_list = beta_list, LP_list = LP_list, eta_list = eta_list)
  return(results)
}
