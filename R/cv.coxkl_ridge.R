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
#' @param lambda Optional numeric scalar or vector of penalty parameters. If `NULL`, a sequence is generated automatically.
#' @param nfolds Integer; number of cross-validation folds. Default = \code{5}.
#' @param cv.criteria Character string specifying the cross-validation criterion. Choices are:
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
                           lambda = NULL, nlambda = 100, lambda.min.ratio = 1e-3, nfolds = 5, 
                           cv.criteria = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"),
                           c_index_stratum = NULL,
                           message = FALSE, seed = NULL,  ...) {
  
  ## Input checks & data preparation
  if (is.null(etas)) stop("etas must be provided.", call. = FALSE)
  etas <- sort(etas)
  cv.criteria <- match.arg(cv.criteria, choices = c("V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"))
  
  if (is.null(RS) && is.null(beta)) {
    stop("No external information is provided. Either RS or beta must be provided.")
  } else if (is.null(RS) && !is.null(beta)) {
    if (length(beta) != ncol(z)) stop("beta dimension mismatch with z.", call. = FALSE)
    RS <- as.matrix(z) %*% as.matrix(beta)
  } else {
    RS <- as.matrix(RS)
  }
  
  if (is.null(stratum)) {
    warning("Stratum not provided. Treating all data as one stratum.", call. = FALSE)
    stratum <- rep(1, nrow(z))
  } else {
    if (!is.null(c_index_stratum) & !identical(stratum, c_index_stratum)) {
      stop("The provided 'c_index_stratum' is not identical to 'stratum'!")
    }
    stratum <- match(stratum, unique(stratum))
  }
  
  ## Sort by (stratum, time)
  time_order <- order(stratum, time)
  time <- as.numeric(time[time_order])
  stratum <- as.numeric(stratum[time_order])
  z <- as.matrix(z)[time_order, , drop = FALSE]
  delta <- as.numeric(delta[time_order])
  RS <- RS[time_order, , drop = FALSE]
  n <- nrow(z)
  
  ## Precompute per-stratum counts for full data (used by VVH/LinPred external)
  n.each_stratum_full <- as.numeric(table(stratum))
  
  ## Fixed CV folds
  if (!is.null(seed)) set.seed(seed)
  folds <- get_fold(nfolds = nfolds, delta = delta, stratum = stratum)
  
  n_eta <- length(etas)
  results_list <- list()
  
  ## External baseline accumulators (criterion-matched)
  if (cv.criteria == "V&VH") {
    pl_full_RS <- pl_cal_theta(as.vector(RS), delta, n.each_stratum_full)
    ext_vvh_per_fold <- numeric(nfolds)
  } else if (cv.criteria == "CIndex_pooled") {
    ext_numer <- numeric(nfolds); ext_denom <- numeric(nfolds)
  } else if (cv.criteria == "CIndex_foldaverage") {
    ext_c_per_fold <- numeric(nfolds)
  }
  ## For LinPred, external baseline = -2 * pl_cal_theta(RS, full data), computed once at end.
  
  if (message) {
    cat("Cross-validation over eta sequence:\n")
    pb_eta <- txtProgressBar(min = 0, max = n_eta, style = 3, width = 30)
  }
  
  ## Outer loop over eta
  for (ei in seq_along(etas)) {
    eta <- etas[ei]
    
    ## Determine lambda sequence for this eta
    if (is.null(lambda)) {
      fit0 <- coxkl_ridge(z = z, delta = delta, time = time, stratum = stratum,
                          RS = RS, eta = eta, lambda = NULL,
                          nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                          data_sorted = TRUE, message = FALSE, ...)
      lambda_seq <- as.vector(fit0$lambda)
    } else {
      lambda_seq <- sort(lambda, decreasing = TRUE)
    }
    L <- length(lambda_seq)
    
    ## Accumulators for this eta over folds
    if (cv.criteria == "V&VH") {
      vvh_sum <- rep(0, L)
    } else if (cv.criteria == "LinPred") {
      Y <- matrix(NA_real_, nrow = n, ncol = L)
      colnames(Y) <- round(lambda_seq, 6)
    } else if (cv.criteria == "CIndex_pooled") {
      numer <- rep(0, L); denom <- rep(0, L)
    } else if (cv.criteria == "CIndex_foldaverage") {
      csum <- rep(0, L); cnt <- rep(0, L)
    }
    
    ## Inner loop over folds
    for (f in seq_len(nfolds)) {
      train_idx <- which(folds != f)
      test_idx  <- which(folds == f)
      
      ## External baseline per fold (criterion-matched) — compute once per fold (independent of eta/lambda)
      if (ei == 1) {  # do this only once across etas; RS is fixed
        if (cv.criteria == "V&VH") {
          n.each_stratum_train <- as.numeric(table(stratum[train_idx]))
          pl_train_RS <- pl_cal_theta(as.vector(RS[train_idx]), delta[train_idx], n.each_stratum_train)
          ext_vvh_per_fold[f] <- pl_full_RS - pl_train_RS
        } else if (cv.criteria == "CIndex_pooled") {
          if (is.null(c_index_stratum)) {
            stratum_test <- stratum[test_idx]
          } else {
            stratum_test <- c_index_stratum[test_idx]
          }
          cstat_ext <- c_stat_stratcox(time[test_idx],
                                       as.vector(RS[test_idx]),
                                       stratum_test,
                                       delta[test_idx])
          ext_numer[f] <- cstat_ext$numer
          ext_denom[f] <- cstat_ext$denom
        } else if (cv.criteria == "CIndex_foldaverage") {
          if (is.null(c_index_stratum)) {
            stratum_test <- stratum[test_idx]
          } else {
            stratum_test <- c_index_stratum[test_idx]
          }
          cstat_ext <- c_stat_stratcox(time[test_idx],
                                       as.vector(RS[test_idx]),
                                       stratum_test,
                                       delta[test_idx])
          ext_c_per_fold[f] <- cstat_ext$c_statistic
        }
      }
      
      ## Fit the entire lambda path on training fold (warm-start handled inside coxkl_ridge)
      fit_f <- coxkl_ridge(z = z[train_idx, , drop = FALSE],
                           delta = delta[train_idx],
                           time = time[train_idx],
                           stratum = stratum[train_idx],
                           RS = RS[train_idx, , drop = FALSE],
                           eta = eta, lambda = lambda_seq,
                           data_sorted = TRUE, message = FALSE, ...)
      
      beta_mat <- fit_f$beta
      z_train <- z[train_idx, , drop = FALSE]
      z_test  <- z[test_idx, , drop = FALSE]
      delta_train <- delta[train_idx]
      stratum_train <- stratum[train_idx]
      
      LP_train <- z_train %*% beta_mat
      LP_all   <- z       %*% beta_mat
      LP_test  <- z_test  %*% beta_mat
      
      ## Fold-wise CV statistics for this eta (vector over lambda grid)
      if (cv.criteria == "V&VH") {
        n.each_stratum_train <- as.numeric(table(stratum_train))
        n.each_stratum_all   <- as.numeric(table(stratum))
        pl_all <- apply(LP_all,  2, function(col) pl_cal_theta(col, delta,        n.each_stratum_all))
        pl_tr  <- apply(LP_train,2, function(col) pl_cal_theta(col, delta_train, n.each_stratum_train))
        vvh_sum <- vvh_sum + (pl_all - pl_tr)
      } else if (cv.criteria == "LinPred") {
        Y[test_idx, ] <- LP_test
      } else {
        if (is.null(c_index_stratum)) {
          stratum_test <- stratum[test_idx]
        } else {
          stratum_test <- c_index_stratum[test_idx]
        }
        ## Column-wise C-index over lambda (no dependency on unavailable c_stat_stratcox_vec)
        numer_vec <- denom_vec <- cstat_vec <- rep(NA_real_, ncol(LP_test))
        for (j in seq_len(ncol(LP_test))) {
          cstat_j <- c_stat_stratcox(time[test_idx], LP_test[, j], stratum_test, delta[test_idx])
          numer_vec[j] <- cstat_j$numer
          denom_vec[j] <- cstat_j$denom
          cstat_vec[j] <- cstat_j$c_statistic
        }
        if (cv.criteria == "CIndex_pooled") {
          numer <- numer + numer_vec
          denom <- denom + denom_vec
        } else {
          csum <- csum + cstat_vec
          cnt  <- cnt  + 1
        }
      }
    } # end folds
    
    ## Aggregate folds for this eta
    if (cv.criteria == "V&VH") {
      cve_eta <- vvh_sum
    } else if (cv.criteria == "LinPred") {
      Lmat <- loss.coxkl_highdim(delta, Y, stratum, total = FALSE)
      cve_eta <- colSums(Lmat)
    } else if (cv.criteria == "CIndex_pooled") {
      cve_eta <- numer / denom
    } else {
      cve_eta <- csum / pmax(cnt, 1)
    }
    
    results_list[[ei]] <- data.frame(
      eta = eta,
      lambda = lambda_seq,
      score = cve_eta,
      stringsAsFactors = FALSE
    )
    
    if (message) setTxtProgressBar(pb_eta, ei)
  } # end etas
  
  if (message) close(pb_eta)
  
  ## Bind per-eta results (full 2D grid)
  results_df <- do.call(rbind, results_list)
  
  ## Label columns by criterion and also compute best-per-eta
  if (cv.criteria %in% c("V&VH", "LinPred")) {
    results_df$Loss <- -2 * results_df$score
    results_df$score <- NULL
    best_per_eta <- do.call(rbind, lapply(split(results_df, results_df$eta), function(df) {
      df[which.min(df$Loss), , drop = FALSE]
    }))
  } else if (cv.criteria == "CIndex_pooled") {
    results_df$CIndex_pooled <- results_df$score
    results_df$score <- NULL
    best_per_eta <- do.call(rbind, lapply(split(results_df, results_df$eta), function(df) {
      df[which.max(df$CIndex_pooled), , drop = FALSE]
    }))
  } else { # CIndex_foldaverage
    results_df$CIndex_foldaverage <- results_df$score
    results_df$score <- NULL
    best_per_eta <- do.call(rbind, lapply(split(results_df, results_df$eta), function(df) {
      df[which.max(df$CIndex_foldaverage), , drop = FALSE]
    }))
  }
  rownames(best_per_eta) <- NULL
  
  ## External baseline (single scalar), matched to the chosen criterion
  external_stat <- switch(cv.criteria,
                          "V&VH" = {
                            # Same aggregation as internal VVH: sum over folds then convert to loss scale
                            -2 * sum(ext_vvh_per_fold)
                          },
                          "LinPred" = {
                            # Assembled CV linear predictor equals RS itself
                            -2 * pl_cal_theta(as.vector(RS), delta, n.each_stratum_full)
                          },
                          "CIndex_pooled" = {
                            sum(ext_numer) / sum(ext_denom)
                          },
                          "CIndex_foldaverage" = {
                            mean(ext_c_per_fold)
                          }
  )
  
  return(list(
    external_stat.full_results = results_df,   # all (eta, lambda) combinations with CV statistic
    external_stat.best_per_eta = best_per_eta, # per-eta best lambda rows
    external_stat = external_stat# scalar baseline matched to cv.criteria
  ))
}

  
  
