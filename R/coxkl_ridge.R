#' Cox Proportional Hazards Model with Ridge Penalty and External Information
#'
#' Fits a Cox proportional hazards model using a ridge-type penalty (L2) on all covariates.
#' The model can integrate external information either as precomputed risk scores (`RS`) 
#' or externally supplied coefficients (`beta`). A tuning parameter `eta` controls the 
#' relative weight of the external information. If `lambda` is not provided, a lambda 
#' sequence is automatically generated.
#'
#' @param z Numeric matrix of covariates (observations in rows, predictors in columns).
#' @param delta Numeric vector of event indicators (1 = event, 0 = censored).
#' @param time Numeric vector of observed times.
#' @param stratum Optional numeric or factor vector specifying strata.
#' @param RS Optional numeric vector or matrix of external risk scores.
#' @param beta Optional numeric vector of externally derived coefficients.
#' @param eta Non-negative scalar controlling the strength of external information.
#' @param lambda Optional numeric scalar or vector of penalty parameters. If `NULL`, a sequence is generated automatically.
#' @param nlambda Number of lambda values to generate if `lambda` is `NULL`.
#' @param lambda.min.ratio Ratio defining the minimum lambda relative to `lambda.max`.
#' @param penalty.factor Numeric scalar in `[0, 1)` .Controls the overall strength of the penalty when generating the ridge 
#'   regression lambda sequence. Smaller values correspond to stronger penalization. Only used when `lambda = NULL`.
#' @param tol Convergence tolerance for the iterative estimation algorithm.
#' @param Mstop Maximum number of iterations for estimation.
#' @param backtrack Logical; if `TRUE`, uses backtracking line search.
#' @param message Logical; if `TRUE`, progress messages are printed during model fitting. Default is `FALSE`.
#' @param data_sorted Logical; if `TRUE`, assumes input data is already sorted by strata and time.
#' @param ... Additional arguments.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{lambda}}{The lambda sequence used for estimation.}
#'   \item{\code{beta}}{Matrix of estimated coefficients for each lambda (columns correspond to lambda).}
#'   \item{\code{linear.predictors}}{Matrix of linear predictors for each lambda.}
#'   \item{\code{likelihood}}{Vector of log-partial likelihoods for each lambda.}
#' }
#'
#' @export
coxkl_ridge <- function(z, delta, time, stratum = NULL, RS = NULL, beta = NULL, eta = NULL,
                        lambda = NULL, nlambda = 100, lambda.min.ratio = 1e-3, penalty.factor = 0.999,
                        tol = 1.0e-4, Mstop = 50, backtrack = FALSE, message = FALSE, data_sorted = FALSE,
                        beta_initial = NULL, ...){
  if (is.null(eta)){
    warning("eta is not provided. Setting eta = 0 (no external information used).", call. = FALSE)
    eta <- 0
  } else {
    if (!is.finite(eta) || eta < 0 || length(eta) != 1) {
      stop("eta must be a non-negative scalar.", call.=FALSE)
    }
  }
  
  if (is.null(RS) && is.null(beta)) {
    stop("Error: No external information is provided. Either RS or beta must be provided.")
  } else if (is.null(RS) && !is.null(beta)) {
    if (length(beta) == ncol(z)) {
      if (message) message("External beta information is used.")
      RS <- as.matrix(z) %*% as.matrix(beta)
    } else {
      stop("Error: The dimension of beta does not match the number of columns in z.")
    }
  } else if (!is.null(RS)) {
    RS <- as.matrix(RS)
    if (message) message("External Risk Score information is used.")
  }
  
  if (!data_sorted) {
    ## ---- Sorting Section ----
    if (is.null(stratum)) {
      warning("Stratum information not provided. All data is assumed to originate from a single stratum!", call. = FALSE)
      stratum <- rep(1, nrow(z))
      time_order <- order(time)
    } else {
      stratum <- match(stratum, unique(stratum))
      time_order <- order(stratum, time)
    }
    
    time <- as.numeric(time[time_order])
    stratum <- as.numeric(stratum[time_order])
    z_mat <- as.matrix(z)[time_order, , drop = FALSE]
    delta <- as.numeric(delta[time_order])
    RS <- as.numeric(RS[time_order, , drop = FALSE])
  } else {
    z_mat <- as.matrix(z)
    time <- as.numeric(time)
    delta <- as.numeric(delta)
    stratum <- as.numeric(stratum)
    RS <- as.numeric(RS)
  }
  
  n.each_stratum <- as.numeric(table(stratum))
  beta.init <- rep(0, ncol(z_mat)) #initial value of beta
  delta_tilde <- calculateDeltaTilde(delta, time, RS, n.each_stratum)
  
  if (is.null(lambda)) {
    if (nlambda < 2) {
      stop("nlambda must be at least 2", call. = FALSE)
    } else if (nlambda != round(nlambda)){
      stop("nlambda must be a positive integer", call. = FALSE)
    }
    # will use lasso lambda sequence
    lambda.fit <- setupLambdaCoxKL(z_mat, time, delta, delta_tilde, RS, beta.init, stratum, 
                                   group = 1:ncol(z_mat), group.multiplier = rep(1, ncol(z_mat)), n.each_stratum, alpha = 1 - penalty.factor,
                                   eta, nlambda, lambda.min.ratio)
    lambda.seq <- lambda.fit$lambda.seq
  } else {
    nlambda <- length(lambda)  # Note: lambda can be a single value
    lambda.seq <- as.vector(sort(lambda, decreasing = TRUE))
  }
  
  LP_mat <- matrix(NA, nrow = nrow(z_mat), ncol = nlambda)
  beta_mat <- matrix(NA, nrow = ncol(z_mat), ncol = nlambda)
  likelihood_mat <- rep(NA, nlambda)
  
  lambda_names <- round(lambda.seq, 4)
  colnames(LP_mat) <- lambda_names
  colnames(beta_mat) <- lambda_names
  names(likelihood_mat) <- lambda_names
  

  
  if (is.null(beta_initial)){
    beta_initial <- rep(0, ncol(z_mat))
  }
  
  if (message) {
    cat("Cross-validation over lambda sequence:\n")
    pb <- txtProgressBar(min = 0, max = nlambda, style = 3, width = 30)
  }
  
  for (i in seq_along(lambda.seq)) {
    lambda <- lambda.seq[i]
    beta_est <- KL_Cox_Estimate_cpp(z_mat, delta, delta_tilde, n.each_stratum, eta, beta_initial,
                                    tol, Mstop, lambda = lambda, backtrack = backtrack, message = FALSE)
    LP_train <- z_mat %*% as.matrix(beta_est)
    beta_mat[, i] <- beta_est
    LP_mat[, i] <- LP_train
    likelihood_mat[i] <- pl_cal_theta(LP_train, delta, n.each_stratum)
    
    beta_initial <- beta_est  # "warm start"
    if (message) setTxtProgressBar(pb, i)
  }
  if (message) close(pb)
  
  
  if (data_sorted == FALSE){
    LinPred_original <- matrix(NA_real_, nrow = length(time_order), ncol = nlambda)
    LinPred_original[time_order, ] <- LP_mat
  } else {
    LinPred_original <- LP_mat
  }
  
  structure(list(
    lambda = lambda.seq,
    beta = beta_mat,
    linear.predictors = LinPred_original,
    likelihood = likelihood_mat),
    class = "coxkl")  
}




