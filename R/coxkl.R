#' Cox Proportional Hazards Model with KL Divergence for Data Integration
#'
#' Fits a Cox proportional hazards model that incorporates external information
#' using Kullback–Leibler (KL) divergence. External information can be supplied
#' either as precomputed external risk scores (`RS`) or as externally derived
#' coefficients (`beta`). The integration strength is controlled by the tuning
#' parameter(s) `eta`.
#'
#' @param z Numeric matrix of covariates with rows representing observations and
#'   columns representing predictor variables. All covariates must be numeric.
#' @param delta Numeric vector of event indicators (1 = event, 0 = censored).
#' @param time Numeric vector of observed event or censoring times. No sorting
#'   required.
#' @param stratum Optional numeric or factor vector defining strata.
#' @param RS Optional numeric vector or matrix of external risk scores. Length
#'   (or number of rows) must equal the number of observations. If not supplied,
#'   `beta` must be provided.
#' @param beta Optional numeric vector of external coefficients (e.g., from prior
#'   studies). Length must equal the number of columns in `z`. Use zeros to
#'   represent covariates without external information. If not supplied, `RS`
#'   must be provided.
#' @param etas Numeric vector of tuning parameters controlling the reliance on
#'   external information. Larger values place more weight on the external
#'   source.
#' @param tol Convergence tolerance for the optimization algorithm. Default is
#'   `1e-4`.
#' @param Mstop Maximum number of iterations for the optimization algorithm.
#'   Default is `100`.
#' @param backtrack Logical; if `TRUE`, backtracking line search is applied during
#'   optimization. Default is `FALSE`.
#' @param message Logical; if `TRUE`, progress messages are printed during model
#'   fitting. Default is `FALSE`.
#' @param data_sorted Logical; if `TRUE`, input data are assumed to be already
#'   sorted by stratum and time. Default is `FALSE`.
#'
#' @return An object of class \code{"coxkl"}, implemented as a list with the
#'   following components:
#'   \describe{
#'     \item{\code{eta}}{Numeric vector of tuning parameters used.}
#'     \item{\code{beta}}{Numeric matrix of estimated coefficients. Columns
#'       correspond to different values of \code{eta}, rows to covariates.}
#'     \item{\code{linear.predictors}}{Numeric matrix of linear predictors on the
#'       training data, aligned with rows of \code{z}. If input data were not
#'       pre-sorted, linear predictors are re-ordered back to the original input
#'       order.}
#'     \item{\code{likelihood}}{Numeric vector of log-partial likelihood values,
#'       one per fitted model (indexed by \code{eta}).}
#'   }
#'
#' @examples
#' data(ExampleData)
#' etas <- generate_eta(method = "exponential", n = 5, max_eta = 5)
#'
#' # Example 1: use external beta information
#' result1 <- coxkl(
#'   z = ExampleData$z,
#'   delta = ExampleData$status,
#'   time = ExampleData$time,
#'   stratum = ExampleData$stratum, 
#'   beta = ExampleData$beta_external,
#'   etas = etas
#' )
#'
#' # Example 2: use external risk score information
#' rs <- as.matrix(ExampleData$z) %*% as.matrix(ExampleData$beta_external)
#' result2 <- coxkl(
#'   z = ExampleData$z,
#'   delta = ExampleData$status,
#'   time = ExampleData$time,
#'   stratum = ExampleData$stratum,
#'   RS = rs,
#'   etas = etas
#' )
#'
#' @export


coxkl <- function(z, delta, time, stratum = NULL,
                  RS = NULL, beta = NULL, 
                  etas, tol = 1.0e-4, Mstop = 100,
                  backtrack = FALSE,
                  message = FALSE,
                  data_sorted = FALSE,
                  beta_initial = NULL){
  
  if (is.null(RS) && is.null(beta)) {
    stop("No external information is provided. Either RS or beta must be provided.")
  } else if (is.null(RS) && !is.null(beta)) {
    if (length(beta) == ncol(z)) {
      if (message) message("External beta information is used.")
      RS <- as.matrix(z) %*% as.matrix(beta)
    } else {
      stop("The dimension of beta does not match the number of columns in z.")
    }
  } else if (!is.null(RS)) {
    RS <- as.matrix(RS)
    if (message) message("External Risk Score information is used.")
  }
  
  if (!data_sorted) {
    ## ---- Sorting Section ----
    if (is.null(stratum)) {
      if (message) warning("Stratum information not provided. All data is assumed to originate from a single stratum!", call. = FALSE)
      stratum <- rep(1, nrow(z))
    } else {
      stratum <- match(stratum, unique(stratum))
    }
    time_order <- order(stratum, time)
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
  
  etas <- sort(etas)
  n_eta <- length(etas)
  LP_mat <- matrix(NA, nrow = nrow(z_mat), ncol = n_eta)
  beta_mat <- matrix(NA, nrow = ncol(z_mat), ncol = n_eta)
  likelihood_mat <- rep(NA, n_eta)
  
  eta_names <- round(etas, 4)
  colnames(LP_mat) <- eta_names
  colnames(beta_mat) <- eta_names
  names(likelihood_mat) <- eta_names
  
  n.each_stratum <- as.numeric(table(stratum))
  delta_tilde <- calculateDeltaTilde(delta, time, RS, n.each_stratum)
  
  if (is.null(beta_initial)){
    beta_initial <- rep(0, ncol(z_mat))
  }
  
  if (message) {
    cat("Cross-validation over eta sequence:\n")
    pb <- txtProgressBar(min = 0, max = n_eta, style = 3, width = 30)
  }
  
  for (i in seq_along(etas)){  #"etas" already in ascending order
    eta <- etas[i]
    beta_est <- KL_Cox_Estimate_cpp(z_mat, delta, delta_tilde, n.each_stratum, eta, beta_initial,
                                    tol, Mstop, lambda = 0, backtrack = backtrack, message = F)
    LP <- z_mat %*% as.matrix(beta_est)
    LP_mat[, i] <- LP
    beta_mat[, i] <- beta_est
    likelihood_mat[i] <- pl_cal_theta(LP, delta, n.each_stratum)
    
    beta_initial <- beta_est  # "warm start"
    if (message) setTxtProgressBar(pb, i)
  }
  if (message) close(pb)
  
  if (data_sorted == FALSE){
    LinPred_original <- matrix(NA_real_, nrow = length(time_order), ncol = n_eta)
    LinPred_original[time_order, ] <- LP_mat
  } else {
    LinPred_original <- LP_mat
  }
  
  structure(list(
    eta = etas,
    beta = beta_mat,
    linear.predictors = LinPred_original,
    likelihood = likelihood_mat),
    class = "coxkl")  
}

