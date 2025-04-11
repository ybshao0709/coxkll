#' Cox Proportional Hazards Model with KL Divergence using Cross-Validation to select best tuning parameter.
#'
#' This function estimates the coefficients of a Cox proportional hazards model using
#' Kullback-Leibler divergence for data integration, allowing for the incorporation of
#' an external risk score (RS) and adjustment of the integration via a list of tuning parameters (eta).
#'
#' @param z A matrix of covariates where rows represent observations and columns represent variables.
#' @param delta A vector of time-to-event outcomes, where 1 indicates an event has occurred and 0 indicates censoring.
#' @param time A vector of event times corresponding to each observation.
#' @param RS A vector or matrix of external risk scores corresponding to each observation, used for data integration.
#' @param eta_list A vector of tuning parameters for data integration, allowing for adjustment in the influence of the external risk scores.
#' @param tol The tolerance level for the convergence criterion of the optimization algorithm. Default is 1.0e-7.
#' @param Mstop The maximum number of iteration steps for the optimization algorithm. Default is 50.
#' @param nfolds The number of cross-validation folds. Default is 5.
#' @param criteria The criterion used to evaluate model performance during
#' cross-validation.  
#' 
#' @return A list containing:
#'   \item{result}{A vector of performance scores (e.g., C-index) for each `eta`
#'     value in `eta_list`. The length of this vector is equal to the length
#'     of `eta_list`.}
#'   \item{eta}{The `eta_list` provided as input.}
#'   \item{criteria}{The criterion used to evaluate model performance during
#' cross-validation.}
#' 
#' @examples
#' data(ExampleData)
#' cv_result <- cv.coxkl(z = z, delta = delta, time = time,
#'                        RS = rs_external1, eta_list = eta_list,
#'                        nfolds = 5, criteria = "C-Index")
#'                        
#' @export
cv.coxkl <- function(z, delta, time, RS, beta, eta_list, tol=1.0e-7, Mstop = 50, nfolds = 5, criteria = "C-Index"){
  
  if(is.null(RS) && is.null(beta)) {
    stop("Error: No external information is provided. Either RS or beta must be provided.")
  }  else if(is.null(RS) && !is.null(beta)) {
    # Check if the dimension of beta matches the number of columns in z
    if(length(beta) == ncol(z)) {
      print("External beta information is used.")
      RS <- as.matrix(Z_internal) %*% as.matrix(beta_external_homo)
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
  
  likelihood_all <- NULL
  
  folds <- get_fold(nfolds = nfolds, delta)
  
  for (eta_index in seq_along(eta_list)){
    
    eta = eta_list[eta_index]
    likelihood_cv = rep(0, nfolds)
    folds <- get_fold(5, delta)
    for(f in 1:nfolds){
      train_idx <- which(folds != f)
      test_idx <- which(folds == f)
      
      # Extract training and testing data
      Z_train <- z[train_idx, ]
      delta_train <- delta[train_idx]
      time_train <- time[train_idx]
      
      Z_test <- z[test_idx, ]
      delta_test <- delta[test_idx]
      time_test <- time[test_idx]
      
      cox_estimate <- coxkl_sorted(z = Z_train, delta = delta_train, time = time_train, 
                                  RS = RS[train_idx], eta_list = eta, tol=tol, Mstop = Mstop)
      
      beta_train <- cox_estimate$beta_list[[1]]
      
      # if (criteria == "V&VH")
      # {
      #   LP_train <- as.matrix(Z_train)%*%as.matrix(beta_train)
      #   LP_internal <- as.matrix(z)%*%as.matrix(beta_train)
      #   
      #   likelihood_cv[cv] <- pl_cal_theta(LP_internal, delta, t_internal) - pl_cal_theta(LP_train, delta_train, t_train)        
      # }
      # 
      #C-Index
      if (criteria == "C-Index"){
        LP_test <- as.matrix(Z_test)%*%as.matrix(beta_train)
        likelihood_cv[f] <- glmnet::Cindex(LP_test, Surv(time_test, delta_test))        
      }
    }
    
    if (criteria == "V&VH"){
      likelihood[k] <- mean(likelihood_cv)
    }
    
    if (criteria == "C-Index"){
      likelihood_all <- c(likelihood_all, mean(likelihood_cv))
    }
    
  }
  results <- list(result = likelihood_all, eta_list = eta_list, criteria = criteria)
  return(results)
}



get_fold <- function(nfolds = 5, delta){
  n <- length(delta)
  ind1 <- which(delta==1)
  ind0 <- which(delta==0)
  n1 <- length(ind1)
  n0 <- length(ind0)
  fold1 <- 1:n1 %% nfolds
  fold0 <- (n1 + 1:n0) %% nfolds
  fold1[fold1==0] <- nfolds
  fold0[fold0==0] <- nfolds
  fold <- integer(n)
  fold[delta==1] <- sample(fold1)
  fold[delta==0] <- sample(fold0)
  return(fold)
}

pl_cal_theta <- function(lp, delta, time){
  delta = delta[order(time)]
  lp = lp[order(time)]
  S0 <- rev(cumsum(rev(exp(lp))))
  pl <- sum(delta*(lp-log(S0)))
  return(pl)
}



# cv.coxkl_lasso <- function(z, delta, time, RS, beta, eta_list, lambda, tol=1.0e-7, Mstop = 50, nfolds = 5, 
#                            penalty="lasso",
#                            alpha=1,
#                            criteria = "V&VH"){
#   
#   if(is.null(RS) && is.null(beta)) {
#     stop("Error: No external information is provided. Either RS or beta must be provided.")
#   }  else if(is.null(RS) && !is.null(beta)) {
#     # Check if the dimension of beta matches the number of columns in z
#     if(length(beta) == ncol(z)) {
#       print("External beta information is used.")
#       RS <- as.matrix(z) %*% as.matrix(beta)
#     } else {
#       stop("Error: The dimension of beta does not match the number of columns in z.")
#     }
#   } else if(!is.null(RS)) {
#     print("External Risk Score information is used.")
#   }
#   
#   time_order <- order(time)
#   delta      <- delta[time_order]
#   z          <- z[time_order,]
#   RS <- RS[time_order]
#   time       <- time[time_order]
#   
#   likelihood_all <- NULL
#   
#   # folds <- get_fold(nfolds = nfolds, delta)
#   
#   # print(folds)
#   
#   for (eta_index in seq_along(eta_list)){
#     
#     eta = eta_list[eta_index]
#     likelihood_cv = rep(0, nfolds)
#     folds <- get_fold(nfolds, delta)
#     for(f in 1:nfolds){
#       train_idx <- sort(which(folds != f))  
#       test_idx <- sort(which(folds == f))   
#       
#       # Extract training and testing data
#       Z_train <- z[train_idx, ]
#       delta_train <- delta[train_idx]
#       time_train <- time[train_idx]
#       
#       Z_test <- z[test_idx, ]
#       delta_test <- delta[test_idx]
#       time_test <- time[test_idx]
#       
#       # cox_estimate <- coxkl_sorted(z = Z_train, delta = delta_train, time = time_train,
#       #                              RS = RS[train_idx], eta_list = eta, tol=tol, Mstop = Mstop)
#       # 
#       # beta_train <- cox_estimate$beta_list[[1]]
#       
#       res_train <- cv.coxkl_highdim(z = Z_train, delta = delta_train, time = time_train,
#                                     theta_tilde = RS[train_idx],
#                                     group = c(1:dim(Z_train)[2]), K = 1,
#                                     penalty = penalty,
#                                     eta_kl = eta,
#                                     alpha = alpha,
#                                     lambda = lambda,
#                                     cv.method = "LinPred")
#       
#       beta_train <- res_train$fit$beta[,which(res_train$lambda == res_train$lambda.min)]
#       
#       if (criteria == "V&VH")
#       {
#         LP_train <- as.matrix(Z_train)%*%as.matrix(beta_train)
#         LP_internal <- as.matrix(z)%*%as.matrix(beta_train)
#         
#         likelihood_cv[f] <- pl_cal_theta(LP_internal, delta_test, time_test) - pl_cal_theta(LP_train, delta_train, time_train)
#       }
#       
#       #C-Index
#       if (criteria == "C-Index"){
#         LP_test <- as.matrix(Z_test)%*%as.matrix(beta_train)
#         likelihood_cv[f] <- glmnet::Cindex(LP_test, Surv(time_test, delta_test))
#       }
#     }
#     
#     if (criteria == "V&VH"){
#       likelihood_all[eta_index] <- mean(likelihood_cv)
#     }
#     
#     if (criteria == "C-Index"){
#       likelihood_all <- c(likelihood_all, mean(likelihood_cv))
#     }
#     
#   }
#   results <- list(result = likelihood_all, eta_list = eta_list, criteria = criteria)
#   return(results)
# }
