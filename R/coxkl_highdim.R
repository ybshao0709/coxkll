#' Cox Proportional Hazards Model with KL Divergence for Data Integration and Lasso & Elastic Net Penalty
#'
#' Fits a Cox proportional hazards model that incorporates external information
#' using Kullback–Leibler (KL) divergence, with an optional L1 (Lasso) or elastic net penalty on
#' the coefficients. External information can be supplied either as precomputed external 
#' risk scores (`RS`) or as externally derived coefficients (`beta`). The integration 
#' strength is controlled by the tuning parameter `eta`.
#'
#' @details
#' Setting `lambda = 0` reduces to the unpenalized \code{coxkl()} model.
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
#' @param eta Numeric tuning parameter controlling the reliance on external
#'   information. Larger values place more weight on the external source.
#' @param alpha Elastic-net mixing parameter in \eqn{(0,1]}. When \eqn{\alpha=1}
#'   the penalty is lasso.
#' @param lambda Optional nonnegative penalty parameter(s). If a numeric vector
#'   is supplied, the path is taken as-is. If `NULL`, a sequence is generated
#'   using `nlambda` and `lambda.min.ratio`.
#' @param nlambda Integer number of lambda values to generate when `lambda` is
#'   `NULL`. Default `100`.
#' @param lambda.min.ratio Ratio of the smallest to the largest lambda when
#'   generating a sequence (when `lambda` is `NULL`). Default `1e-3`.
#' @param lambda.early.stop Logical; if `TRUE`, stop traversing the lambda path
#'   early based on convergence or screening criteria. Default `FALSE`.
#' @param tol Convergence tolerance for the optimization algorithm. Default is
#'   `1e-3`.
#' @param Mstop Maximum number of iterations for the inner optimization at a
#'   given lambda. Default is `1000`.
#' @param max.total.iter Maximum total iterations across the entire lambda path.
#'   Default is `(Mstop * nlambda)`.
#' @param group Integer vector of group indices defining group
#'   membership of predictors for grouped penalties; use `0` to indicate
#'   unpenalized variables.
#' @param group.multiplier A vector of values representing multiplicative factors 
#'   by which each covariate's penalty is to be multiplied. Default is a vector of 1's.
#' @param standardize Logical; if `TRUE`, columns of `z` are standardized prior
#'   to fitting, with coefficients re-scaled on output. Default `TRUE`.
#' @param nvar.max Integer cap on the number of active variables allowed during
#'   fitting. Default number of predictors.
#' @param group.max Integer cap on the number of active groups allowed during
#'   fitting. Default total number of groups.
#' @param stop.loss.ratio Relative improvement threshold for early stopping along
#'   the path; optimization may stop if objective gain falls below this value.
#'   Default `1e-3`.
#' @param actSet Logical; if `TRUE`, use an active-set strategy. Default `TRUE`.
#' @param actIter Maximum number of active-set refinement iterations per lambda.
#'   Default `Mstop`.
#' @param actGroupNum Maximum number of active groups allowed under the
#'   active-set scheme. 
#' @param actSetRemove Logical; if `TRUE`, allow dropping variables/groups from
#'   the active set during iterations. Default `FALSE`.
#' @param trace.lambda Logical; if `TRUE`, record path-wise traces across the
#'   lambda sequence. Default `FALSE`.
#' @param message Logical; if `TRUE`, progress messages are printed during model
#'   fitting. Default is `FALSE`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class \code{"coxkl_highdim"}, a list with components:
#' \describe{
#'   \item{\code{beta}}{Coefficient estimates (vector or matrix across the path).}
#'   \item{\code{group}}{A \code{factor} of the original group assignments.}
#'   \item{\code{lambda}}{The lambda value(s) used or generated.}
#'   \item{\code{alpha}}{The elastic-net mixing parameter used.}
#'   \item{\code{likelihood}}{Vector of log-partial likelihoods for each lambda.}
#'   \item{\code{n}}{Number of observations.}
#'   \item{\code{df}}{Effective degrees of freedom (e.g., number of nonzero
#'     coefficients or group-adjusted count) along the path.}
#'   \item{\code{iter}}{Number of iterations taken (per lambda and/or total).}
#'   \item{\code{W}}{Exponentiated linear predictors on the original scale}.}
#'   \item{\code{group.multiplier}}{Group-specific penalty multipliers used.}
#' }
#'
#' @examples
#' data(ExampleData)
#'
#' # With external beta
#' fit_hd1 <- coxkl_highdim(
#'   z = ExampleData$z,
#'   delta = ExampleData$status,
#'   time = ExampleData$time,
#'   stratum = ExampleData$stratum,
#'   beta = ExampleData$beta_external,
#'   etas = 1
#' )
#'
#'
#' @export
coxkl_highdim <- function(z, delta, time, stratum = NULL, RS = NULL, beta = NULL, eta = NULL,
                          alpha = NULL, lambda = NULL, nlambda = 100, lambda.min.ratio = 1e-3, lambda.early.stop = FALSE,
                          tol = 1.0e-4, Mstop = 1000, max.total.iter = (Mstop * nlambda), 
                          group = 1:ncol(z), group.multiplier = NULL, standardize = T, 
                          nvar.max = ncol(z), group.max = length(unique(group)), stop.loss.ratio = 1e-3, 
                          actSet = TRUE, actIter = Mstop, actGroupNum = sum(unique(group) != 0), actSetRemove = F,
                          returnX = FALSE, trace.lambda = FALSE, message = FALSE, data_sorted = FALSE, ...){
  
  if (is.null(alpha)){
    warning("alpha is not provided. Setting alpha = 1 (lasso penalty).", call. = FALSE)
    alpha <- 1
  } else if (alpha > 1 | alpha <= 0) {
    stop("alpha must be in (0, 1]", call.=FALSE)
  }
  
  if (is.null(eta)){
    warning("eta is not provided. Setting eta = 0 (no external information used).", call. = FALSE)
    eta <- 0
  } else {
    if (!is.finite(eta) || eta < 0 || length(eta) != 1) {
      stop("eta must be a non-negative scalar.", call.=FALSE)
    }
  }
  
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
  
  z <- as.matrix(z)
  delta <- as.numeric(delta)
  time <- as.numeric(time)
  
  
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
    z <- as.matrix(z)[time_order, , drop = FALSE]
    delta <- as.numeric(delta[time_order])
    RS <- as.numeric(RS[time_order, , drop = FALSE])
  } else {
    z <- as.matrix(z)
    time <- as.numeric(time)
    delta <- as.numeric(delta)
    stratum <- as.numeric(stratum)
    RS <- as.numeric(RS)
  }
  
  n.each_stratum <- as.numeric(table(stratum))
  
  initial.group <- group
  if (standardize == T){
    std.Z <- newZG.Std(z, group, group.multiplier)
  } else {
    std.Z <- newZG.Unstd(z, group, group.multiplier)
  }
  Z <- std.Z$std.Z[, , drop = F]
  group <- std.Z$g  
  group.multiplier <- std.Z$m 
  
  p <- ncol(Z)
  n <- length(delta)
  nvar.max <- as.integer(nvar.max)
  group.max <- as.integer(group.max)

  beta.init <- rep(0, ncol(Z)) #initial value of beta
  delta_tilde <- calculateDeltaTilde(delta, time, RS, n.each_stratum)
  
  if (is.null(lambda)) {
    if (nlambda < 2) {
      stop("nlambda must be at least 2", call. = FALSE)
    } else if (nlambda != round(nlambda)){
      stop("nlambda must be a positive integer", call. = FALSE)
    }
    lambda.fit <- setupLambdaCoxKL(Z, time, delta, delta_tilde, RS, beta.init, stratum, 
                                   group, group.multiplier, n.each_stratum, alpha,
                                   eta, nlambda, lambda.min.ratio)
    lambda.seq <- lambda.fit$lambda.seq
    beta <- lambda.fit$beta
  } else {
    nlambda <- length(lambda)  # Note: lambda can be a single value
    lambda.seq <- as.vector(sort(lambda, decreasing = TRUE))
    beta <- beta.init
  }
  
  K <- as.integer(table(group)) #number of features in each group
  K0 <- as.integer(if (min(group) == 0) K[1] else 0)
  K1 <- as.integer(if (min(group) == 0) cumsum(K) else c(0, cumsum(K)))
  
  initial.active.group <- -1
  if (actSet == TRUE){
    if (K0 == 0){
      initial.active.group <- which(K == min(K))[1] - 1
    }
  } else {
    actIter <- Mstop
  }
  
  fit <- KL_Cox_highdim(Z, delta, delta_tilde, eta, n.each_stratum, beta, K1, K0, 
                        lambda.seq, alpha, lambda.early.stop, stop.loss.ratio, 
                        group.multiplier, max.total.iter, Mstop, tol, 
                        initial.active.group, nvar.max, group.max, trace.lambda, 
                        actSet, actIter, actGroupNum, actSetRemove)
  # colSums(fit$beta != 0)   #internal check for non-zero coefficients (when at lambda_max, beta should be all zeros)
  
  beta <- fit$beta
  LinPred <- fit$LinPred
  df <- fit$Df
  iter <- fit$iter
  loss <- fit$loss
  
  # Eliminate saturated lambda values
  ind <- !is.na(iter)
  lambda <- lambda.seq[ind]
  beta <- beta[, ind, drop = FALSE]
  loss <- loss[ind]
  LinPred <- LinPred[, ind, drop = FALSE]
  df <- df[ind]
  iter <- iter[ind]
  
  if (iter[1] == max.total.iter){
    stop("Algorithm failed to converge for any values of lambda", call. = FALSE)
  }
  if (sum(iter) == max.total.iter){
    warning("Algorithm failed to converge for all values of lambda", call. = FALSE)
  }
  
  
  # Original scale
  beta <- unorthogonalize(beta, std.Z$std.Z, group)
  rownames(beta) <- colnames(Z)
  if (std.Z$reorder == TRUE){  # original order of beta
    beta <- beta[std.Z$ord.inv, , drop = F]
  }
  if (standardize == T) {
    original.beta <- matrix(0, nrow = length(std.Z$scale), ncol = ncol(beta))
    original.beta[std.Z$nz, ] <- beta / std.Z$scale[std.Z$nz]
    beta <- original.beta
  }

  
  # Names
  dimnames(beta) <- list(colnames(Z), round(lambda, digits = 4))
  colnames(LinPred) <- round(lambda, digits = 4)
  
  #recover the original order of linear predictors
  if (data_sorted == FALSE){
    LinPred_original <- matrix(NA_real_, nrow = length(time_order), ncol = ncol(LinPred))
    LinPred_original[time_order, ] <- LinPred
  } else {
    LinPred_original <- LinPred
  }
  

  result <- structure(list(beta = beta,
                           group = factor(initial.group),
                           lambda = lambda,
                           alpha = alpha,
                           likelihood = loss,
                           n = n,
                           df = df,
                           iter = iter,
                           W = exp(LinPred_original),  # rescale beta will not change the linear predictors (Z matrix is also standardized)
                           df = df,
                           group.multiplier = group.multiplier),
                      class = "coxkl_highdim")
  
  if (returnX == TRUE){  # used for cross validation!
    result$returnX <- list(XX = std.Z,
                           time = time,
                           delta = delta,
                           stratum = stratum,
                           RS = RS)
  }
  return(result)
}
