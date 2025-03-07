
cv.migrlasso <- function(X, y, K, ..., nfolds=10, seed, fold, se=c('quick', 'bootstrap'), 
                            cv.method = c('LinPred', 'VVH'),
                            returnY=FALSE, trace=FALSE) {
  
  se <- match.arg(se)
  # Complete data fit
  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$K <- K
  fit.args$returnX <- TRUE
  #fit <- do.call("grpsurv", fit.args)
  fit <- do.call("migrlasso", fit.args)
  # Get standardized X, y
  X <- fit$XG
  y <- cbind(fit$time, fit$fail)
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX) fit$X <- NULL
  
  # Set up folds
  n <- nrow(X)
  if (!missing(seed)) set.seed(seed)
  if (missing(fold)) {
    ind1 <- which(fit$fail==1)
    ind0 <- which(fit$fail==0)
    n1 <- length(ind1)
    n0 <- length(ind0)
    fold1 <- 1:n1 %% nfolds
    fold0 <- (n1 + 1:n0) %% nfolds
    fold1[fold1==0] <- nfolds
    fold0[fold0==0] <- nfolds
    fold <- integer(n)
    fold[fit$fail==1] <- sample(fold1)
    fold[fit$fail==0] <- sample(fold0)
  } else {
    nfolds <- max(fold)
  }
  Y <- matrix(NA, nrow=n*K, ncol=length(fit$lambda))

  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$warn <- FALSE
  
  if(cv.method == 'LinPred'){
    for (i in 1:nfolds) {
      #if (trace) cat("Starting CV fold #", i, sep="","\n")
      res <- cvf.surv(i, X, y, K, fold, cv.args)
      index.foldi <- which(fold==i)
      index.foldi2 <- NULL
      for (j in 1:K) {
        index.foldi2 <- c(index.foldi2, (index.foldi+n*(j-1)))
      }
      Y[index.foldi2, 1:res$nl] <- res$yhat
    }
    # Eliminate saturated lambda values, if any
    # ind <- which(apply(is.finite(Y), 2, all))
    # print(ind)
    # Y <- Y[, ind]
    # lambda <- fit$lambda[ind]
    # Return
    #if (se == "quick") {
      L <- loss.MIGroupLasso(y, Y, K, total=FALSE)
      cve <- apply(L, 2, sum)/sum(fit$fail)
      #cvse <- apply(L, 2, sd)*sqrt(nrow(L))/sum(fit$fail)
    #} else {
      #cve <- as.double(loss.grpsurv(y, Y))/sum(fit$fail)
      #cvse <- se.grpsurv(y, Y)/sum(fit$fail)
    #}
  } else {
    lambda  <- cv.args$lambda
    cve     <- rep(0, length(lambda))
    for (i in 1:nfolds) {
      if (trace) cat("Starting CV fold #", i, sep="","\n")
      res                   <- cvf.survVVH(i, X, y, K, fold, cv.args)
      EtaAll                <- res$yhatall
      EtaK                  <- res$yhat
      L_K                   <- loss.MIGroupLasso(y[fold!=i,], EtaK, K, total=FALSE)
      L_all                 <- loss.MIGroupLasso(y, EtaAll, K, total=FALSE)
      cve                   <- cve +  (apply(L_all, 2, sum)/sum(fit$fail) - apply(L_K, 2, sum)/sum(fit$fail))
    }
  }
  
  min <- which.min(as.double(cve))
  
  val <- list(cve=cve, fold=fold, lambda=cv.args$lambda, fit=fit, min=min, lambda.min=cv.args$lambda[min], null.dev=cve[1])
  if (returnY) val$Y <- Y
  structure(val, class=c("cv.grpsurv", "cv.grpreg"))
}

cvf.surv <- function(i, XX, y, K, fold, cv.args) {
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i,]
  cv.args$K <- K
  fit.i <- do.call("migrlasso", cv.args)
  
  X2 <- XX[fold==i, , drop=FALSE]
  y2 <- y[fold==i,]
  nl <- length(fit.i$lambda)
  yhat <- predict.MIGroupLasso(fit.i, X2, K)
  
  list(nl=length(fit.i$lambda), yhat=yhat)#, loss=loss)
}


cvf.survVVH <- function(i, XX, y, K, fold, cv.args) {
  cv.args$X <- XX[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i,]
  cv.args$K <- K
  fit.i     <- do.call("migrlasso", cv.args)
  #
  X2        <- XX[fold!=i, , drop=FALSE]
  y2        <- y[fold!=i,]
  nl        <- length(fit.i$lambda)
  yhat      <- predict.MIGroupLasso(fit.i, X2, K)
  yhatall    <- predict.MIGroupLasso(fit.i, XX, K)
  
  list(nl=length(fit.i$lambda), yhat=yhat, yhatall=yhatall)#, loss=loss)
}


