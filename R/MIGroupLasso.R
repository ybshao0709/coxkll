migrlasso <- function(X, y, K, 
                         alpha=1, nlambda=100, lambda,
                         lambda.min={if (nrow(X) > ncol(X)) 0.001 else .05}, eps=.001, max.iter=10000,
                         dfmax=p*K, gmax=p, tau=1/3,
                         warn=TRUE, returnX=FALSE, ...) {
  
  #check
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha > 1 | alpha <= 0) stop("alpha must be in (0, 1]", call.=FALSE)
  
  gamma = 3
  
  # Construct XG, Y
  Y <- newS(y)
  group <- NULL
  for (i in 1:K) {
    group <- c(group, c(1:p))
  }
  std <- standardize(X)
  XX <- std[[1]]
  center <- as.vector(std[[2]])
  scale <- std[[3]]
  XX <- orthogonalize(XX, group)
  
  if (nrow(XX) != length(Y$fail)) stop("X and y do not have the same number of observations", call.=FALSE)
  
  
  # Set up lambda
  if (missing(lambda)) {
    lambda <- setupLambda_MIGrLasso(XX, Y$time, Y$fail, K, alpha, lambda.min, nlambda)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  
  #fit
  n <- length(Y$time)
  p <- ncol(XX)/K
  # K <- as.integer(table(XG$g))
  # K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)                     # =0 if each variable has a group index
  # K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))  
  dfmax <- p*K
  warn  <- FALSE
  res   <- MIGrLasso_cox(XX, Y$fail, p = p, K = K, lambda = lambda, 
                         alpha = alpha, eps = eps, max_iter = max.iter, dfmax = dfmax, gmax = p, 
                         warn = warn, user = user.lambda)  
  
  b     <- matrix((res$beta),p*K,nlambda)
  iter  <- as.vector(res$iter)
  df    <- as.vector(res$df)
  loss  <- -1*as.vector(res$Loss)
  Eta   <- matrix(t(res$Eta),n,nlambda)
  j_selected <- res$j_selected
  
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  df <- df[ind]
  loss <- loss[ind]
  #if (iter[1] == max.iter) stop("Algorithm failed to converge for any values of lambda.  This indicates a combination of (a) an ill-conditioned feature matrix X and (b) insufficient penalization.  You must fix one or the other for your model to be identifiable.", call.=FALSE)
  if (warn & any(iter==max.iter)) warning("Algorithm failed to converge for all values of lambda", call.=FALSE)
  
  
  # Unstandardize
  b <- unorthogonalize(b,XX, group, intercept=FALSE)
  beta <- matrix(0,length(scale), ncol=ncol(b))
  nz <- which(scale > 1e-6)
  beta[nz,] <- b / scale[nz]
  xnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
  dimnames(beta) <- list(xnames, round(lambda, digits=4))
  
  # Output
  val <- structure(list(beta = beta,
                        group = group,
                        lambda = lambda,
                        gamma = gamma,
                        alpha = alpha,
                        #loss = loss,
                        #n = n,
                        df = df,
                        iter = iter,
                        time = Y$time,
                        fail = Y$fail,
                        W = exp(Eta),
                        class = c("grpsurv", "grpreg"),
                        j_selected = j_selected))
  if (returnX) {
    val$XG <- XX
  }
  val
}






