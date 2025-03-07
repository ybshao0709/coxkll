predict.MIGroupLasso <- function(object, X, K, 
                            type=c("link", "response", "survival", "median", "norm", 
                                   "coefficients", "vars", "nvars", "groups", "ngroups"),
                            lambda, which=1:length(object$lambda), ...) {
  type <- match.arg(type)

  if (!missing(lambda)) {
    ind <- approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    x <- ind %% 1
    beta <- (1-x)*object$beta[, l, drop=FALSE] + x*object$beta[, r, drop=FALSE]
    colnames(beta) <- round(lambda, 4)
  } else {
    beta <- object$beta[, which, drop=FALSE]
  }
  
  N   <- dim(X)[1]
  p   <- dim(X)[2]/K
  eta <- matrix(0,N*K,length(object$lambda))
  for (i in 1:K) {
    eta[(N*(i-1)+1) : (N*i),] <- X[, (p*(i-1)+1) : (p*i)] %*% beta[(p*(i-1)+1) : (p*i),]
    #eta[(N*(i-1)+1) : (N*i),] <- X[, (p*(1-1)+1) : (p*1)] %*% beta[(p*(1-1)+1) : (p*1),]
  }
  
  return(drop(eta))
}
