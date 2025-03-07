
setupLambda_MIGrLasso <- function(X, y, Delta, K, alpha, lambda.min, nlambda) {
  n <- nrow(X)

  
  ## Fit to unpenalized covariates
  # K <- table(group)
  # K1 <- as.integer(if (min(group)==0) cumsum(K) else c(0, cumsum(K)))
  # if (K1[1]!=0) {
  #   SURV <- get("Surv", asNamespace("survival"))
  #   COXPH <- get("coxph", asNamespace("survival"))
  #   nullFit <- COXPH(SURV(y, Delta) ~ X[, group==0, drop=FALSE])
  #   eta <- nullFit$linear.predictors
  #   rsk <- rev(cumsum(rev(exp(eta))))
  #   s <- Delta - exp(eta)*cumsum(Delta/rsk)
  # } else {
  w <- 1/(n-(1:n)+1)
  s <- Delta - cumsum(Delta*w)
  s <- rep(s, K)
  # }
    
  ## Determine lambda.max
  zmax <- maxgrad(X,s,K,p)
  
  lambda.max <- zmax/alpha
  
  if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)), 0)
  else lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len=nlambda))
  #lambda*sqrt(K)
}

