#' Cross-Validation for CoxKL Ridge Model (eta tuning)
#'
#' This function performs cross-validation on the Cox model with Kullback–Leibler (KL) 
#' penalty and ridge (L2) regularization. It tunes the parameter \code{eta} 
#' (external information weight) using user-specified cross-validation criteria, 
#' while internally selecting the optimal \code{lambda} at each fold.
#'
#' @param z Numeric matrix of covariates with rows representing individuals and
#'   columns representing predictors.
#' @param delta Numeric vector of event indicators (1 = event, 0 = censored).
#' @param time Numeric vector of observed times (event or censoring).
#' @param stratum Optional factor or numeric vector indicating strata.
#' @param RS Optional numeric vector or matrix of external risk scores. If not provided,
#'   \code{beta} must be supplied.
#' @param beta Optional numeric vector of external coefficients (length equal to
#'   \code{ncol(z)}). If not provided, \code{RS} must be supplied.
#' @param etas Numeric vector of candidate \code{eta} values to be evaluated.
#' @param nfolds Integer; number of cross-validation folds. Default = \code{5}.
#' @param cv.eta.criteria Character string specifying the cross-validation criterion
#'   for selecting \code{eta}. Choices are:
#'   \itemize{
#'     \item \code{"V&VH"} (default): V&VH loss.
#'     \item \code{"LinPred"}: loss based on cross-validated linear predictors approach.
#'     \item \code{"CIndex_pooled"}: pool all held-out predictions and compute one overall CIndex.
#'     \item \code{"CIndex_foldaverage"}: average CIndex across folds.
#'   }.
#'   Default is \code{"V&VH"}.
#' @param c_index_stratum Optional stratum vector. Required only when 
#'   \code{cv.eta.criteria} is set to \code{"CIndex_pooled"} or \code{"CIndex_foldaverage"}, 
#'   and a stratified C-index needs to be computed while the fitted model 
#'   is non-stratified. Default is \code{NULL}.
#' @param cv.lambda.criteria Character string specifying the criterion for selecting optimal
#'   \code{lambda}. Choices are:
#'   \itemize{
#'     \item \code{"V&VH"} (default): V&VH loss.
#'     \item \code{"LinPred"}: loss based on cross-validated linear predictors approach.
#'   }
#' @param message Logical; whether to print progress messages. Default = \code{FALSE}.
#' @param seed Optional integer random seed for fold assignment.
#' @param ... Additional arguments passed to \code{\link{coxkl_ridge}}.
#'
#' @return A \code{data.frame} with two columns:
#'   \describe{
#'     \item{\code{eta}}{Candidate \code{eta} values.}
#'     \item{criteria column}{The value of the chosen \code{cv.eta.criteria}
#'       for each \code{eta}. Column name depends on the criterion:
#'       \code{"VVH_Loss"}, \code{"LinPred_Loss"}, \code{"CIndex_pooled"},
#'       or \code{"CIndex_foldaverage"}.}
#'   }
#'   
#' @export


cv.coxkl_ridge <- function(z, delta, time, stratum = NULL, RS = NULL, beta = NULL, etas,
                           nfolds = 5, cv.eta.criteria = "V&VH", c_index_stratum = NULL,
                           cv.lambda.criteria = "V&VH",
                           message = FALSE, seed = NULL,  ...) {
  
  if (is.null(etas)){
    stop("etas must be provided.", call. = FALSE)
  }
  
  cv.eta.criteria <- match.arg(cv.eta.criteria, choices = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"))
  cv.lambda.criteria <- match.arg(cv.lambda.criteria, choices = c("V&VH", "LinPred"))
  
  if (cv.eta.criteria != cv.lambda.criteria){
    warning("The Cross-Validation criteria for eta and lambda are different!")
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
    
    
    if (cv.eta.criteria == "V&VH") {
      cv_vvh <- numeric(nfolds)
    } else if (cv.eta.criteria == "LinPred") {
      cv_all_linpred <- matrix(NA, nrow = length(delta), ncol = 1)
    } else if (cv.eta.criteria == "CIndex_pooled") {
      cv_pooled_cindex_mat <- matrix(0, nfolds, 2)
    } else if (cv.eta.criteria == "CIndex_foldaverage") {
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
      cv.lambda.args <- list(...)  #other input arguments for "coxkl_highdim"
      # cv.lambda.args <- list()
      cv.lambda.args$z <- z_train
      cv.lambda.args$delta <- delta_train
      cv.lambda.args$time <- time_train
      cv.lambda.args$stratum <- stratum_train
      cv.lambda.args$RS <- RS_train
      cv.lambda.args$eta <- eta  #current eta
      cv.lambda.args$nfolds <- nfolds  
      cv.lambda.args$seed <- seed
      cv.lambda.args$cv.lambda.criteria <- cv.lambda.criteria
      cv.lambda.args$trace.cv <- FALSE
      
      invisible(
        capture.output(
          fit <- do.call("cv.lambda.ridge", cv.lambda.args)
        )
      )
      
      beta_train <- fit$best.beta
      
      if (cv.eta.criteria == "V&VH") {
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
        if (cv.eta.criteria == "LinPred") {
          cv_all_linpred[test_idx,] <- LP_test
        } else { #C-Index
          if (is.null(c_index_stratum)){
            stratum_test <- stratum[test_idx]
          } else {
            stratum_test <- c_index_stratum[test_idx]
          }
          if (cv.eta.criteria == "CIndex_pooled") {
            cstat <- c_stat_stratcox(time_test, LP_test, stratum_test, delta_test)
            cv_pooled_cindex_mat[f,] <- c(cstat$numer, cstat$denom)
          } else if (cv.eta.criteria == "CIndex_foldaverage") {
            cstat <- c_stat_stratcox(time_test, LP_test, stratum_test, delta_test)$c_statistic
            cv_cindex[f] <- cstat
          }
        }
      }
    }
    
    if (cv.eta.criteria == "V&VH"){
      result_vec[eta_index] <- sum(cv_vvh)
    } else if (cv.eta.criteria == "LinPred"){
      result_vec[eta_index] <- pl_cal_theta(cv_all_linpred, delta, as.numeric(table(stratum)))
    } else if (cv.eta.criteria == "CIndex_foldaverage") {
      result_vec[eta_index] <- mean(cv_cindex)
    } else if (cv.eta.criteria == "CIndex_pooled"){
      result_vec[eta_index] <- sum(cv_pooled_cindex_mat[, 1], na.rm = TRUE) / sum(cv_pooled_cindex_mat[, 2], na.rm = TRUE)
    }
  }
  
  results <- data.frame(eta = etas)
  if (cv.eta.criteria == "V&VH") {
    results$VVH_Loss <- -2 * result_vec
  } else if (cv.eta.criteria == "LinPred") {
    results$LinPred_Loss <- -2 * result_vec
  } else if (cv.eta.criteria == "CIndex_pooled") {
    results$CIndex_pooled <- result_vec
  } else if (cv.eta.criteria == "CIndex_foldaverage") {
    results$CIndex_foldaverage <- result_vec
  }
  return(results)
}


