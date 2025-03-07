#' Cox Proportional Hazards Model with KL Divergence for High-Dimensional Data
#'
#' This function estimates the coefficients of a Cox proportional hazards model
#' using Kullback-Leibler (KL) divergence for data integration in high-dimensional
#' settings. It incorporates an external risk score (represented by `theta_tilde`)
#' and supports regularization via lasso or ridge penalties.  It is designed for
#' situations where the number of covariates (p) can be larger than the number
#' of observations (n).
#'
#' @param z A matrix of covariates. Rows represent observations, and columns
#'   represent variables.
#' @param delta A vector of event indicators (1 = event, 0 = censored).
#' @param time A vector of survival/censoring times.
#' @param group A vector indicating group membership for each covariate.  Used
#'   for group penalties (not implemented in this documentation, but present in
#'   the function signature).
#' @param K The number of groups.
#' @param penalty The type of penalty to apply.  Must be either "lasso" (for
#'   L1 regularization) or "ridge" (for L2 regularization).
#' @param theta_tilde A vector representing the external risk score information.
#'    This is *crucially* different from the `RS` parameter in the other `coxkl`
#'    functions. Here, `theta_tilde` is used directly in the KL divergence
#'    calculation, *not* as a pre-calculated risk score.  It's likely related
#'    to coefficients from an external model.
#' @param eta_kl A *single* numeric value controlling the influence of the
#'   external risk score (`theta_tilde`). A higher `eta_kl` gives more weight
#'   to the external information. Default is 0.
#' @param alpha The elastic net mixing parameter.  `alpha = 1` corresponds to
#'   lasso, `alpha = 0` corresponds to ridge, and values between 0 and 1
#'   represent a mixture of the two.  Default is 1 (lasso).
#' @param nlambda The number of lambda values to use in the regularization path.
#'   Default is 100.  The function automatically generates a sequence of lambda
#'   values.
#' @param lambda An optional user-specified sequence of lambda values.  If
#'   provided, `nlambda` is ignored.  It's generally recommended to let the
#'   function generate the lambda sequence unless you have specific requirements.
#' @param lambda.min The ratio of the smallest lambda value to the largest lambda
#'   value in the automatically generated sequence.  Default is 0.001.
#' @param eps Convergence threshold for the coordinate descent algorithm.
#'   Default is 0.001.
#' @param max.iter The maximum number of iterations for the coordinate descent
#'   algorithm. Default is 1000.
#' @param dfmax Maximum number of non-zero coefficients.
#' @param gmax Maximum number of groups with non-zero coefficients (relevant
#'   for group penalties, but not fully described here).
#' @param tau A parameter related to the penalty (specific meaning depends on
#'   the penalty type, but not fully described here). Default is 1/3.
#' @param group.multiplier An optional vector of weights for each group, used
#'  to apply different penalties to different groups (not fully described here).
#' @param warn Logical value indicating whether to print warnings. Default is
#'   TRUE.
#' @param returnX Logical value. If TRUE, the design matrix `z` is included
#' in the returned list.
#' @param activeSet Logical. If TRUE, an active set strategy is used to speed up computation.
#' @param actIter Integer. Maximum number of iteration in active set cycling.
#' @param actNum Integer. Maximum number of variables added to the active set at each step.
#' @param ... Additional arguments (not used in this simplified documentation,
#'   but potentially used by the underlying implementation).
#'
#' @return A list containing:
#'    \item{beta}{A matrix of coefficients.  Each *column* corresponds to a
#'      different value of `lambda` in the regularization path.  The number of
#'      rows is equal to the number of covariates (number of columns in `z`).}
#'    \item{lambda}{The sequence of `lambda` values used.}
#'    \item{df}{The number of non-zero coefficients for each value of `lambda`.}
#'    \item{iter}{The number of iterations taken for each value of `lambda`.}
#'    \item{group}{The `group` vector provided as input.}
#'   \item{penalty}{The penalty used}
#'   \item{alpha}{The `alpha` provided as input.}
#'    \item{X}{If `returnX = TRUE`, the design matrix `z` is returned.}
#'
#' @export
#'
#' @examples
#' # Load example data (assuming the package is installed and loaded)
#' data(simulatedData)
#' result_lasso <- coxkl_highdim(z = z, delta = delta, time = time, group = group,
#'                              K = K, penalty = "lasso", theta_tilde = theta_tilde,
#'                              eta_kl = 1)
#'                              
#'                            
coxkl_highdim <-function(z, delta, time, group, K, penalty="lasso",
                         theta_tilde,
                         eta_kl = 0, 
                         alpha=1, nlambda=100, lambda,
                         lambda.min=0.001, eps=.001, max.iter=1000,
                         dfmax=p, gmax=length(unique(group)), tau=1/3, group.multiplier,
                         warn=TRUE, returnX=FALSE, 
                         activeSet = FALSE,
                         actIter=50,
                         actNum=5,
                         ...){
  gamma = 3
  p = ncol(z)
  n = length(delta)
  bilevel <- FALSE
  
  z <- as.matrix(z)
  delta_matrix <- as.matrix(delta)
  time_matrix <- as.matrix(time)
  
  std <- standardize(z)
  XX <- std[[1]]
  center <- as.vector(std[[2]])
  scale <- std[[3]]
  
  K <- as.integer(table(group))
  group.multiplier <- rep(rep(sqrt(K),2),p)
  
  if (missing(lambda)) {
    lambda <- setupLambdaCox_self(XX, time_matrix, delta_matrix, group, penalty, alpha, lambda.min, nlambda, group.multiplier)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  
  # group <- c(1)
  K <- as.integer(table(group))
  K0 <- as.integer(if (min(group)==0) K[1] else 0)          # =0 if each variable has a group index
  K1 <- as.integer(if (min(group)==0) cumsum(K) else c(0, cumsum(K))) 
  dfmax <- p
  warn <- TRUE
  
  delta_tilde <- calculateDeltaTilde(delta_matrix, time_matrix, theta_tilde)
  
  res <- gdfit_cox_kl(X = z, d = delta_matrix, penalty = penalty,
                      delta_tilde = delta_tilde,
                      eta_kl = eta_kl,
                      K1 = K1, K0 = K0,
                      lambda = lambda,
                      alpha = alpha, eps = eps, max_iter = max.iter,
                      gamma = as.double(gamma), group_multiplier = group.multiplier,
                      dfmax = as.integer(dfmax),
                      gmax = length(unique(group)),
                      warn = as.integer(warn), user = FALSE)
  
  b     <- matrix(t(res$beta),p*K,nlambda)
  iter  <- as.vector(res$iter)
  df    <- as.vector(res$df)
  loss  <- -1*as.vector(res$Loss)
  Eta   <- matrix(t(res$Eta),n*p,nlambda)
  
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  df <- df[ind]
  loss <- loss[ind]
  #if (iter[1] == max.iter) stop("Algorithm failed to converge for any values of lambda.  This indicates a combination of (a) an ill-conditioned feature matrix X and (b) insufficient penalization.  You must fix one or the other for your model to be identifiable.", call.=FALSE)
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda", call.=FALSE)
  
  
  # Unstandardize
  #print(group)
  #print(dim(b))
  #print(dim(XX))
  #b <- unorthogonalize(b, XX, group, intercept=FALSE)
  #if (XG$reorder) b <- b[XG$ord.inv,]
  #beta <- matrix(0, nrow=length(scale), ncol=ncol(b))
  #beta[XG$nz,] <- b / XG$scale[XG$nz]
  
  #dimnames(beta) <- list(XG$names, round(lambda, digits=4))
  
  # Output
  val <- structure(list(beta = b,
                        group = group,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        loss = loss,
                        n = n,
                        df = df,
                        iter = iter,
                        group.multiplier = group.multiplier,
                        time = time,
                        delta = delta,
                        #fail = Y$fail,
                        W = exp(Eta)),
                   class = c("grpsurv", "grpreg"))
  if (returnX) {
    val$XX <- XX
    val$delta_matrix <- delta_matrix
    val$time_matrix <- time_matrix
  }
  val
}



#' Cox Proportional Hazards Model with KL Divergence
#'
#' This function estimates the coefficients of a Cox proportional hazards model using
#' Kullback-Leibler divergence for data integration, allowing for the incorporation of
#' an external risk score (RS) and adjustment of the integration via a list of tuning parameters (eta).
#' @export
cv.coxkl_highdim <-function(z, delta, time, group,
                            theta_tilde,
                            eta_kl = 0, ..., nfolds = 10, seed, fold, se=c('quick', 'bootstrap'),
                            cv.method = c('LinPred', 'VVH'),
                            returnY=FALSE, trace=FALSE){
  
  se <- match.arg(se)
  p <- ncol(z)
  K <- 1
  
  fit.args <- list(...)
  fit.args$z <- z
  fit.args$delta <- delta
  fit.args$time <- time
  fit.args$group <- c(1:p)
  fit.args$K <- 1
  fit.args$theta_tilde <- theta_tilde
  fit.args$eta_kl <- eta_kl
  
  fit.args$returnX <- TRUE
  fit <- do.call("coxkl_highdim", fit.args)
  
  X <- fit$XX
  delta_matrix <- fit$delta_matrix
  time_matrix <- fit$time_matrix
  returnX <- list()$returnX
  if (is.null(returnX) || !returnX) fit$X <- NULL
  
  # Set up folds
  n <- nrow(X)
  #if (!missing(seed)) set.seed(seed)
  #if (missing(fold)) {
  
  ind1 <- which(fit$delta==1)
  ind0 <- which(fit$delta==0)
  n1 <- length(ind1)
  n0 <- length(ind0)
  fold1 <- 1:n1 %% nfolds
  fold0 <- (n1 + 1:n0) %% nfolds
  fold1[fold1==0] <- nfolds
  fold0[fold0==0] <- nfolds
  fold_sub <- integer(n/K)
  fold <- integer(n) 
  fold_sub[fit$delta==1] <- sample(fold1)
  fold_sub[fit$delta==0] <- sample(fold0)
  fold <- rep(fold_sub, each = K)
  #} else {
  #  nfolds <- max(fold)
  #}
  Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))
  
  cv.args <- list()
  cv.args$lambda <- fit$lambda
  cv.args$group <- fit$group
  #cv.args$group.multiplier <- fit$
  cv.args$warn <- FALSE
  cv.args$eta_kl <- eta_kl
  
  if(cv.method == 'LinPred'){
    for (i in 1:nfolds) {
      res <- cvf.surv(i, z, theta_tilde, fit$delta, fit$time, K, fold, fold_sub, cv.args)
      Y[fold==i, 1:res$nl] <- res$yhat
    }
    ind <- which(apply(is.finite(Y), 2, all))
    lambda <- fit$lambda[ind]
    L <- loss.grpsurv_self_multi(delta_matrix, Y, total=FALSE)
    cve <- apply(L, 2, sum)/sum(delta_matrix)
    # cvse <- apply(L, 2, sd)*sqrt(nrow(L))/sum(delta_matrix)
  } else if (cv.method == 'VVH'){
    lambda  <- cv.args$lambda
    cve     <- rep(0, length(lambda))
    for (i in 1:nfolds) {
      res                   <- cvf.survVVH(i, z, theta_tilde, fit$delta, fit$time, K, fold, fold_sub, cv.args)
      EtaAll                <- res$yhatall
      EtaK                  <- res$yhat
      L_K                   <- loss.grpsurv_self_multi(delta_matrix[fold!=i], EtaK, K, total=FALSE)
      L_all                 <- loss.grpsurv_self_multi(fit$delta, EtaAll, K, total=FALSE)
      cve                   <- cve +  (apply(L_all, 2, sum)/sum(delta_matrix) - apply(L_K, 2, sum)/sum(delta_matrix))
    }
  }
  
  min <- which.min(cve)
  
  
  val <- list(cve=cve, fold=fold, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min], null.dev=cve[1])
  if (returnY) val$Y <- Y
  structure(val, class=c("cv.grpsurv", "cv.grpreg"))
  
}



cvf.surv <- function(i, z, theta_tilde, delta, time, K, fold, fold_sub, cv.args) {
  K = 1
  Xlist_temp <- z[fold_sub!=i, , drop=FALSE]
  
  cv.args$z <- Xlist_temp
  cv.args$time <- time[fold_sub!=i]
  cv.args$delta <- delta[fold_sub!=i]
  cv.args$theta_tilde <- theta_tilde[fold_sub!=i]
  cv.args$K <- K
  fit.i <- do.call("coxkl_highdim", cv.args)
  
  Xlist_temp2 <- z[fold_sub==i, , drop=FALSE]
  
  nl <- length(fit.i$lambda)
  yhat <- predict.grpsurv_self_multi(fit.i, Xlist_temp2, K)
  
  list(nl=length(fit.i$lambda), yhat=yhat)#, loss=loss)
}


cvf.survVVH <- function(i, z, theta_tilde, delta, time, K, fold, fold_sub, cv.args) {
  Xlist_temp <- NULL
  
  Xlist_temp <- z[fold_sub!=i, , drop=FALSE]
  
  cv.args$z <- Xlist_temp
  cv.args$time <- time[fold_sub!=i]
  cv.args$delta <- delta[fold_sub!=i]
  cv.args$theta_tilde <- theta_tilde[fold_sub!=i]
  cv.args$K <- K
  fit.i <- do.call("coxkl_highdim", cv.args)
  
  Xlist_temp2 <- z[fold_sub==i, , drop=FALSE]
  
  nl        <- length(fit.i$lambda)
  yhat      <- predict.grpsurv_self_multi(fit.i, Xlist_temp2, K)
  yhatall   <- predict.grpsurv_self_multi(fit.i, z, K)
  
  list(nl=length(fit.i$lambda), yhat=yhat, yhatall=yhatall)#, loss=loss)
}
