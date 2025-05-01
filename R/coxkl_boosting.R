#' #' Cox Proportional Hazards Model with KL Divergence (Boosting)
#' #'
#' #' This function estimates the coefficients of a Cox proportional hazards model using
#' #' Kullback-Leibler divergence for data integration, allowing for the incorporation of
#' #' an external risk score (RS) and adjustment of the integration via a list of tuning parameters (eta).
#' #'
#' #' @param z A matrix of covariates where rows represent observations and columns represent variables.
#' #' @param delta A vector of time-to-event outcomes, where 1 indicates an event has occurred and 0 indicates censoring.
#' #' @param time A vector of event times corresponding to each observation.
#' #' @param RS A vector or matrix of external risk scores corresponding to each observation, used for data integration.
#' #' @param eta_list A vector of tuning parameters for data integration, allowing for adjustment in the influence of the external risk scores.
#' #' @param tol The tolerance level for the convergence criterion of the optimization algorithm. Default is 1.0e-7.
#' #' @param Mstop The maximum number of iteration steps for the optimization algorithm. Default is 50.
#' #' @param rate The learning rate (step size) for the boosting algorithm.
#' #'   Default is 0.001.  Smaller values typically lead to slower but potentially
#' #'   more stable convergence.
#' #' @param beta 
#' #'   
#' #' @return A list containing three elements:
#' #'   - `betas`: A list of estimated coefficient matrices, one for each value of eta.
#' #'   - `LPs_train`: A list of linear predictors for the training set, one for each value of eta.
#' #'   - `eta`: The original list of eta values provided as input.
#' #'   
#' #' @examples
#' #' # Load example data (assuming the package is installed and loaded)
#' #' data(simulatedData)
#' #' # Fit the model with a single eta value
#' #' result1 <- coxkl_boosting(z = z, delta = delta, time = time,
#' #'                  RS = rs_external1, eta_list = c(0,1,2,3,4))
#' #' 
#' #' @export
#' coxkl_boosting <- function(z, delta, time, RS, eta_list, tol=1e-6, Mstop = 10000, rate = 0.001, beta = NULL){
#'   
#'   if(is.null(RS) && is.null(beta)) {
#'     stop("Error: No external information is provided. Either RS or beta must be provided.")
#'   }  else if(is.null(RS) && !is.null(beta)) {
#'     # Check if the dimension of beta matches the number of columns in z
#'     if(length(beta) == ncol(z)) {
#'       print("External beta information is used.")
#'       RS <- as.matrix(Z_internal) %*% as.matrix(beta_external_homo)
#'     } else {
#'       stop("Error: The dimension of beta does not match the number of columns in z.")
#'     }
#'   } else if(!is.null(RS)) {
#'     print("External Risk Score information is used.")
#'   }
#'   
#'   time_order <- order(time)
#'   delta      <- delta[time_order]
#'   z          <- z[time_order,]
#'   RS <- RS[time_order]
#'   time       <- time[time_order]
#'   
#'   z_mat <- as.matrix(z)
#'   delta_mat <- as.matrix(delta)
#'   p     <- ncol(z_mat)
#'   beta  <- as.matrix(rep(0,p))
#'   RS    <- as.matrix(RS)
#'   
#'   beta_list <- list()
#'   LP_list <- list()
#'   
#'   for (i in seq_along(eta_list)){
#'     eta <- eta_list[i]
#'     res <- klcox_boosting(z_mat, delta, RS, eta = eta,
#'                           rate = rate, tol = tol, maxit = Mstop)
#'     beta_train <- res$beta
#'     LP_train <- as.matrix(z) %*% as.matrix(beta_train)
#' 
#'     beta_list[[i]] <- beta_train
#'     LP_list[[i]] <- LP_train
#'   }
#'   
#'   results <- list(beta_list = beta_list, LP_list = LP_list, eta_list = eta_list)
#'   return(results)
#' }
