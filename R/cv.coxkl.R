#' Cross-Validated Cox–KL to Tune the Integration Parameter (eta)
#'
#' Performs K-fold cross-validation to select the integration parameter `eta`
#' for the Cox–KL model. Each fold fits the model on a training split and
#' evaluates on the held-out split using the specified performance criterion.
#'
#' @param z Numeric matrix of covariates (rows = observations, columns = variables).
#' @param delta Numeric vector of event indicators (1 = event, 0 = censored).
#' @param time Numeric vector of observed event or censoring times.
#' @param stratum Optional numeric or factor vector defining strata.
#' @param RS Optional numeric vector or matrix of external risk scores. If omitted,
#'   `beta` must be supplied.
#' @param beta Optional numeric vector of external coefficients. If omitted, `RS`
#'   must be supplied.
#' @param etas Numeric vector of candidate tuning values to be cross-validated.
#'   Default is `NULL`, which sets `etas = 0`.
#' @param tol Convergence tolerance for the optimizer used inside `coxkl`. Default `1e-4`.
#' @param Mstop Maximum number of Newton iterations used inside `coxkl`. Default `100`.
#' @param backtrack Logical; if `TRUE`, backtracking line search is applied during
#'   optimization. Default is `FALSE`.
#' @param nfolds Number of cross-validation folds. Default `5`.
#' @param criteria Character string specifying the performance criterion.
#'   Choices are `"V&VH"`, `"LinPred"`, `"CIndex_pooled"`, or `"CIndex_foldaverage"`.
#'   Default `"V&VH"`.
#' @param c_index_stratum Optional stratum vector. Only required when
#'   \code{criteria} is set to `"CIndex_pooled"` or `"CIndex_foldaverage"`,
#'   and a stratified C-index is desired while the fitted model is non-stratified.
#'   Default `NULL`.
#' @param message Logical; if `TRUE`, prints progress messages. Default `FALSE`.
#' @param seed Optional integer seed for reproducible fold assignment. Default `NULL`.
#' @param ... Additional arguments passed to \code{\link{coxkl}}.
#'
#' @return A \code{data.frame} with one row per candidate `eta` and columns:
#' \describe{
#'   \item{\code{eta}}{The candidate `eta` values.}
#'   \item{\code{VVH_Loss}}{If \code{criteria = "V&VH"}, the cross-validated V&VH loss.}
#'   \item{\code{LinPred_Loss}}{If \code{criteria = "LinPred"}, the loss based on linear predictors.}
#'   \item{\code{CIndex_pooled}}{If \code{criteria = "CIndex_pooled"}, the pooled cross-validated C-index.}
#'   \item{\code{CIndex_foldaverage}}{If \code{criteria = "CIndex_foldaverage"}, the average fold-wise C-index.}
#' }
#'
#' @examples
#' data(ExampleData)
#' etas <- generate_eta(method = "exponential", n = 5, max_eta = 5)
#'
#' # Example 1: use external beta information, evaluate by V&VH
#' result1 <- cv.coxkl(z = ExampleData$z, delta = ExampleData$status,
#'                     time = ExampleData$time, beta = ExampleData$beta_external,
#'                     etas = etas,
#'                     nfolds = 5, criteria = "V&VH")
#'
#' # Example 2: use external risk score information, evaluate by pooled C-index
#' rs <- as.matrix(ExampleData$z) %*% as.matrix(ExampleData$beta_external)
#' result2 <- cv.coxkl(z = ExampleData$z, delta = ExampleData$status,
#'                     time = ExampleData$time, RS = rs,
#'                     etas = etas,
#'                     nfolds = 5, criteria = "CIndex_pooled")
#'
#' @export

cv.coxkl <- function(z, delta, time, stratum = NULL,
                     RS = NULL, beta = NULL,
                     etas = NULL,
                     tol = 1.0e-4, Mstop = 100,
                     backtrack = FALSE,
                     nfolds = 5,
                     criteria = "V&VH",
                     c_index_stratum = NULL,
                     message = FALSE,
                     seed = NULL,...) {
  criteria <- match.arg(criteria, choices = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"))
  
  if (is.null(etas)){
    stop("etas must be provided.", call. = FALSE)
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
  
  if (is.null(stratum)) {
    warning("Stratum information not provided. All data is assumed to originate from a single stratum!", 
            call. = FALSE)
    stratum <- rep(1, nrow(z))
  } else {
    if (!is.null(c_index_stratum) & !identical(stratum, c_index_stratum)) {
      stop("The provided 'c_index_stratum' is not identical to 'stratum'!")
    }
    stratum <- match(stratum, unique(stratum))
  }
  
  # data is sorted for all & train & test data
  time_order <- order(stratum, time)
  time <- as.numeric(time[time_order])
  stratum <- as.numeric(stratum[time_order])
  z <- as.matrix(z)[time_order, , drop = FALSE]
  delta <- as.numeric(delta[time_order])
  RS <- RS[time_order, , drop = FALSE]
  
  
  n <- nrow(z)
  n_eta <- length(etas)
  result_vec <- numeric(n_eta)
  
  for (eta_index in seq_along(etas)){
    eta <- etas[eta_index]
    if (message) message(sprintf("(%d/%d) eta = %g starts...", eta_index, n_eta, round(eta, 5)))
    
    
    if (criteria == "V&VH") {
      cv_vvh <- numeric(nfolds)
    } else if (criteria == "LinPred") {
      cv_all_linpred <- matrix(NA, nrow = length(delta), ncol = 1)
    } else if (criteria == "CIndex_pooled") {
      cv_pooled_cindex_mat <- matrix(0, nfolds, 2)
    } else if (criteria == "CIndex_foldaverage") {
      cv_cindex <- numeric(nfolds)
    }
    
    if (!is.null(seed)) {
      set.seed(seed)
    } else {
      set.seed(NULL)
    }
    
    if (!is.null(c_index_stratum)){
      folds <- get_fold(nfolds = nfolds, delta, c_index_stratum) #always use stratum information for train-test spliting
    } else {
      folds <- get_fold(nfolds = nfolds, delta, stratum)
    }
    
    for (f in seq_len(nfolds)) {
      if (message) message(sprintf("  CV %d/%d", f, nfolds))
      
      train_idx <- which(folds != f)
      test_idx  <- which(folds == f)
      
      # Training data
      z_train <- z[train_idx, , drop = FALSE]
      delta_train <- delta[train_idx]
      time_train <- time[train_idx]
      stratum_train <- stratum[train_idx]
      RS_train <- RS[train_idx, , drop = FALSE]
      
      # Fit model on training data
      cox_estimate <- coxkl(z = z_train,
                            delta = delta_train,
                            time = time_train,
                            stratum = stratum_train,
                            RS = RS_train,
                            etas = eta,
                            tol = tol,
                            Mstop = Mstop,
                            backtrack = backtrack,
                            message = FALSE, 
                            data_sorted = TRUE)
      
      beta_train <- cox_estimate$beta  #cv needs beta-hat
      
      if (criteria == "V&VH") {
        LP_train <- as.matrix(z_train) %*% as.matrix(beta_train)
        LP_internal <- as.matrix(z) %*% as.matrix(beta_train)
        
        n.each_stratum_train <- as.numeric(table(stratum_train))
        n.each_stratum_internal <- as.numeric(table(stratum))
        
        cv_vvh[f] <- pl_cal_theta(LP_internal, delta, n.each_stratum_internal) - 
          pl_cal_theta(LP_train, delta_train, n.each_stratum_train)
      } else {
        z_test <- z[test_idx, , drop = FALSE]
        delta_test <- delta[test_idx]
        time_test <- time[test_idx]
        LP_test <- as.matrix(z_test) %*% as.matrix(beta_train)
        if (criteria == "LinPred") {
          cv_all_linpred[test_idx,] <- LP_test
        } else { #C-Index
          if (is.null(c_index_stratum)){
            stratum_test <- stratum[test_idx]
          } else {
            stratum_test <- c_index_stratum[test_idx]
          }
          if (criteria == "CIndex_pooled") {
            cstat <- c_stat_stratcox(time_test, LP_test, stratum_test, delta_test)
            cv_pooled_cindex_mat[f,] <- c(cstat$numer, cstat$denom)
          } else if (criteria == "CIndex_foldaverage") {
            cstat <- c_stat_stratcox(time_test, LP_test, stratum_test, delta_test)$c_statistic
            cv_cindex[f] <- cstat
          }
        }
      }
    }
    
    if (criteria == "V&VH"){
      result_vec[eta_index] <- sum(cv_vvh)
    } else if (criteria == "LinPred"){
      result_vec[eta_index] <- pl_cal_theta(cv_all_linpred, delta, as.numeric(table(stratum)))
    } else if (criteria == "CIndex_foldaverage") {
      result_vec[eta_index] <- mean(cv_cindex)
    } else if (criteria == "CIndex_pooled"){
      result_vec[eta_index] <- sum(cv_pooled_cindex_mat[, 1], na.rm = TRUE) / sum(cv_pooled_cindex_mat[, 2], na.rm = TRUE)
    }
  }
  
  results <- data.frame(eta = etas)
  if (criteria == "V&VH") {
    results$VVH_Loss <- -2 * result_vec
  } else if (criteria == "LinPred") {
    results$LinPred_Loss <- -2 * result_vec
  } else if (criteria == "CIndex_pooled") {
    results$CIndex_pooled <- result_vec
  } else if (criteria == "CIndex_foldaverage") {
    results$CIndex_foldaverage <- result_vec
  }
  return(results)
}