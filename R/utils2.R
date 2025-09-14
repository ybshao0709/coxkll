get_fold <- function(nfolds = 5, delta, stratum) {
  n <- length(delta)
  fold <- integer(n)
  
  stratum <- as.factor(stratum)
  strata_levels <- levels(stratum)
  
  for (s in strata_levels) {
    idx <- which(stratum == s)
    delta_s <- delta[idx]
    
    ind1 <- which(delta_s == 1)
    ind0 <- which(delta_s == 0)
    
    n1 <- length(ind1)
    n0 <- length(ind0)
    
    fold1 <- 1:n1 %% nfolds
    fold0 <- (n1 + 1:n0) %% nfolds
    fold1[fold1 == 0] <- nfolds
    fold0[fold0 == 0] <- nfolds
    
    fold_s <- integer(length(idx))
    fold_s[ind1] <- sample(fold1)
    fold_s[ind0] <- sample(fold0)
    
    fold[idx] <- fold_s
  }
  
  return(fold)
}


#' Generate a Sequence of Tuning Parameters (eta)
#'
#' Produces a numeric vector of `eta` values to be used in Cox–KL model.
#'
#' @param method Character string specifying the method to generate `eta`.
#'   Options are `"linear"` for a linearly spaced sequence, or `"exponential"`
#'   for an exponentially spaced sequence scaled to `max_eta`. Default is `"exponential"`.
#' @param n Integer, the number of `eta` values to generate. Default is 10.
#' @param max_eta Numeric, the maximum value of `eta` in the sequence. Default is 5.
#'
#' @return Numeric vector of length `n` containing the generated `eta` values.
#'
#' @examples
#' # Generate 10 exponentially spaced eta values up to 5
#' generate_eta(method = "exponential", n = 10, max_eta = 5)
#'
#' # Generate 5 linearly spaced eta values up to 3
#' generate_eta(method = "linear", n = 5, max_eta = 3)
#'
#' @export
generate_eta <- function(method = "exponential", n = 10, max_eta = 5) {
  if (method == "linear") {
    eta_values <- seq(min_eta, max_eta, length.out = n)
  } else if (method == "exponential") {
    eta_values <- exp(seq(log(1), log(100), length.out = n))
    eta_values <- (eta_values - min(eta_values)) / (max(eta_values) - min(eta_values)) * max_eta
  }
  return(eta_values)
}



c_stat_stratcox <- function(time, xbeta, stratum, delta) {
  stratum <- factor(stratum)
  
  ord <- order(stratum, time)
  time <- time[ord]
  xbeta <- xbeta[ord]
  stratum <- stratum[ord]
  delta <- delta[ord]
  
  # Only evaluate strata with at least 2 individuals
  stratum_sizes <- table(stratum)
  valid_strata <- names(stratum_sizes[stratum_sizes > 1])
  
  # Compute per-stratum concordance
  count <- sapply(valid_strata, function(s) {
    idx <- stratum == s
    cox_c_index(time[idx], xbeta[idx], delta[idx])
  })
  
  numer <- Reduce(`+`, count["numer", ])
  denom <- Reduce(`+`, count["denom", ])
  c_stat <- numer / denom
  
  return(list(numer = numer, 
              denom = denom, 
              c_statistic = c_stat))
}


#' Setup Lambda Sequence for Cox–KL Model (Internal)
#'
#' Generates a sequence of penalty parameters (`lambda`) and initializes coefficients
#' for the Cox–KL model. The maximum lambda (`lambda.max`) is defined as the 
#' smallest value that shrinks all coefficients (`beta`) to zero. Note that different 
#' values of `eta` will produce different `lambda` sequences.
setupLambdaCoxKL <- function(Z, time, delta, delta_tilde, RS, beta.init, stratum,
                             group, group.multiplier, n.each_stratum, alpha, 
                             eta, nlambda, lambda.min.ratio) {
  n <- nrow(Z)
  K <- table(group)
  K1 <- as.integer(if (min(group)==0) cumsum(K) else c(0, cumsum(K)))
  storage.mode(K1) <- "integer"
  if (K1[1]!=0) { ## some covariates are not penalized
    nullFit <- coxkl(Z[, group == 0, drop = FALSE], delta, time, stratum, RS, beta = NULL, eta)
    LinPred <- nullFit$linear.predictors[[1]]
    beta.init <- c(nullFit$beta[[1]], rep(0, length(beta.init) - length(nullFit$beta[[1]])))
    rsk <- c()  
    for (i in seq_along(unique(stratum))){
      rsk <- c(rsk, rev(cumsum(rev(exp(LinPred[stratum == i])))))
    }
    r <- (delta + eta * delta_tilde)/(1 + eta) - exp(LinPred) * cumsum(delta / rsk)
  } else { ## all covariates are penalized
    w <- c()
    h <- c()
    for (i in seq_along(unique(stratum))){
      temp.w <- 1 / (n.each_stratum[i] - (1:n.each_stratum[i]) + 1)
      w <- c(w, temp.w)
      h <- c(h, cumsum(delta[stratum == i] * temp.w))
    }
    r <- (delta + eta * delta_tilde)/(1 + eta) - h
    beta.init <- beta.init
  }
  
  ## Determine lambda.max
  zmax <- maxgrad(Z, r, K1, as.double(group.multiplier)) / n
  lambda.max <- zmax/alpha
  
  if (lambda.min.ratio == 0){
    lambda <- c(exp(seq(log(lambda.max), log(.001 * lambda.max), len = nlambda-1)), 0)
  } else {
    lambda <- exp(seq(log(lambda.max), log(lambda.min.ratio * lambda.max), len = nlambda))
  }
  lambda[1] <- lambda[1] + 1e-5
  ls <- list(beta = beta.init, lambda.seq = lambda)
  return(ls)
}



#---------------------------- Functions for high-dim LASSO ----------------------------#
# set up group information
# m: group multiplier, default (NULL) is the square root of group size of remaining features
setupG <- function(group, m){
  group.factor <- factor(group)
  if (any(levels(group.factor) == '0')) {
    g <- as.integer(group.factor) - 1
    lev <- levels(group.factor)[levels(group.factor) != '0']
  } else {
    g <- as.integer(group.factor)
    lev <- levels(group.factor)
  }
  if (is.numeric(group) | is.integer(group)) {
    lev <- paste0("G", lev)
  }
  if (is.null(m)) {
    m <- rep(NA, length(lev))
    names(m) <- lev
  } else {
    TRY <- try(as.integer(group) == g)
    if (inherits(TRY, 'try-error') || any(!TRY)) stop('Attempting to set group.multiplier is ambiguous if group is not a factor', call. = FALSE)
    if (length(m) != length(lev)) stop("Length of group.multiplier must equal number of penalized groups", call. = FALSE)
    if (storage.mode(m) != "double") storage.mode(m) <- "double"
    if (any(m < 0)) stop('group.multiplier cannot be negative', call.=FALSE)
  }
  # "g" contains the group index of each column, but convert "character" group name into integer
  structure(g, levels = lev, m = m)
}

# remove constant columns if necessary
subsetG <- function(g, nz) { # nz: index of non-constant features
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  new <- g[nz] # only include non-constant columns
  dropped <- setdiff(g, new)  # If the entire group has been dropped
  if (length(dropped) > 0) {
    lev <- lev[-dropped] # remaining group
    m <- m[-dropped]
    group.factor <- factor(new) #remaining group factor
    new <- as.integer(group.factor) - 1 * any(levels(group.factor) == '0')  #new group index
  }
  structure(new, levels = lev, m = m)
}

# reorder group index of features
reorderG <- function(g, m) {
  og <- g
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  if (any(g == 0)) {
    g <- as.integer(relevel(factor(g), "0")) - 1
  }
  if (any(order(g) != 1:length(g))) {
    reorder <- TRUE
    gf <- factor(g)
    if (any(levels(gf) == "0")) {
      gf <- relevel(gf, "0")
      g <- as.integer(gf) - 1
    } else {
      g <- as.integer(gf)
    }
    ord <- order(g)
    ord.inv <- match(1:length(g), ord)
    g <- g[ord]
  } else {
    reorder <- FALSE
    ord <- ord.inv <- NULL
  }
  structure(g, levels = lev, m = m, ord = ord, ord.inv = ord.inv, reorder = reorder)
}

#Feather-level standardization
standardize.Z <- function(Z){
  mysd <- function(z){
    sqrt(sum((z - mean(z))^2)/length(z))
  }
  new.Z <- scale(as.matrix(Z), scale = apply(as.matrix(Z), 2, mysd))
  center.Z <- attributes(new.Z)$`scaled:center`
  scale.Z <- attributes(new.Z)$`scaled:scale`
  new.Z <- new.Z[, , drop = F]
  res <- list(new.Z = new.Z, center.Z = center.Z, scale.Z = scale.Z)
  return(res)
}

## converting standardized betas back to original variables
unstandardize <- function(beta, gamma, std.Z){
  original.beta <- matrix(0, nrow = length(std.Z$scale), ncol = ncol(beta))
  original.beta[std.Z$nz, ] <- beta / std.Z$scale[std.Z$nz]  # modified beta
  original.gamma <- t(apply(gamma, 1, function(x) x - crossprod(std.Z$center, original.beta))) # modified intercepts (gamma)
  return(list(gamma = original.gamma, beta = original.beta))
}


# Group-level orthogonalization (column in new order, from group_0 to group_max)
orthogonalize <- function(Z, group) {
  z.names <- colnames(Z)
  n <- nrow(Z)
  J <- max(group)
  QL <- vector("list", J)
  orthog.Z <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  colnames(orthog.Z) <- z.names
  # unpenalized group will not be orthogonalized
  orthog.Z[, which(group == 0)] <- Z[, which(group == 0)]
  
  # SVD and generate orthogonalized X
  for (j in seq_along(integer(J))) {
    ind <- which(group == j)
    if (length(ind) == 0) { # skip 0-length group
      next
    }
    SVD <- svd(Z[, ind, drop = FALSE], nu = 0)  # Q matrix (orthonormal matrix of eigenvectors)
    r <- which(SVD$d > 1e-10)  #remove extremely small singular values
    QL[[j]] <- sweep(SVD$v[, r, drop = FALSE], 2, sqrt(n)/SVD$d[r], "*") # Q * Lambda^{-1/2}
    orthog.Z[, ind[r]] <- Z[, ind] %*% QL[[j]]  # group orthogonalized X, where (X^T * X)/n = I
  }
  nz <- !apply(orthog.Z == 0, 2, all)  #find all zero
  orthog.Z <- orthog.Z[, nz, drop = FALSE]
  attr(orthog.Z, "QL") <- QL
  attr(orthog.Z, "group") <- group[nz]
  return(orthog.Z)
}

# convert orthogonalized beta back to original scales
unorthogonalize <- function(beta, Z, group) {
  ind <- !sapply(attr(Z, "QL"), is.null)
  QL <- Matrix::bdiag(attr(Z, "QL")[ind]) #block diagonal matrix
  if (sum(group == 0) > 0){ #some groups are unpenalized
    ind0 <- which(group==0)
    original.beta <- as.matrix(rbind(beta[ind0, , drop = FALSE], QL %*% beta[-ind0, , drop = FALSE]))
  } else {  # all groups are penalized
    original.beta <- as.matrix(QL %*% beta)
  }
  return(original.beta)
}



# standardize + orthogonalize covariate matrix
newZG.Std <- function(Z, g, m){
  if (any(is.na(Z))){
    stop("Missing data (NA's) detected in covariate matrix!", call. = FALSE)
  }
  if (length(g) != ncol(Z)) {
    stop ("Dimensions of group is not compatible with Z", call. = FALSE)
  }
  G <- setupG(g, m) # setup group
  std <- standardize.Z(Z)
  std.Z <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  
  small_scales <- which(scale <= 1e-6)
  if (length(small_scales) > 0) {
    stop(
      paste0(
        "The following variables have (near) constant columns: ",
        paste(names(scale)[small_scales], collapse = ", ")
      )
    )
  }
  
  nz <- which(scale > 1e-6)   # non-constant columns
  if (length(nz) != ncol(Z)) {
    std.Z <- std.Z[, nz, drop = F]
    G <- subsetG(G, nz)
  }
  # Reorder groups
  G <- reorderG(G, attr(G, 'm'))
  if (attr(G, 'reorder')){
    std.Z <- std.Z[, attr(G, 'ord')]
  }
  # Group-level orthogonalization
  std.Z <- orthogonalize(std.Z, G)
  g <- attr(std.Z, "group")
  # Set group multiplier if missing
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- sqrt(table(g[g != 0]))
  }
  res <- list(std.Z = std.Z, g = g, m = m, reorder = attr(G, 'reorder'), nz = nz,
              ord.inv = attr(G, 'ord.inv'), center = center, scale = scale)
  return(res)
}

# Only orthogonalize covariate matrix
newZG.Unstd <- function(Z, g, m){
  if (any(is.na(Z))){
    stop("Missing data (NA's) detected in covariate matrix!", call. = FALSE)
  }
  if (length(g) != ncol(Z)) {
    stop ("Dimensions of group is not compatible with Z", call. = FALSE)
  }
  G <- setupG(g, m)
  mysd <- function(x){
    sqrt(sum((x - mean(x))^2)/length(x))
  }
  scale <- apply(as.matrix(Z), 2, mysd)
  
  
  small_scales <- which(scale <= 1e-6)
  if (length(small_scales) > 0) {
    stop(
      paste0(
        "The following variables have (near) constant columns: ",
        paste(names(scale)[small_scales], collapse = ", ")
      )
    )
  }
  
  nz <- which(scale > 1e-6) #remove constant columns
  if (length(nz) != ncol(Z)) {
    std.Z <- Z[, nz, drop = F]
    G <- subsetG(G, nz)
  } else {
    std.Z <- Z
  }
  G <- reorderG(G, attr(G, 'm'))
  if (attr(G, 'reorder')){
    std.Z <- std.Z[, attr(G, 'ord')]
  }
  std.Z <- orthogonalize(std.Z, G)
  g <- attr(std.Z, "group")
  # Set group multiplier if missing
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- sqrt(table(g[g != 0]))
  }
  res <- list(std.Z = std.Z, g = g, m = m, reorder = attr(G, 'reorder'),
              ord.inv = attr(G, 'ord.inv'), nz = nz)
  return(res)
}






#---------------------------- Functions for cv.coxkl_highdim ----------------------------#
coef.coxkl_highdim <- function(fit, lambda = NULL, which = 1:length(fit$lambda), drop = TRUE, ...) {
  if (!is.null(lambda)) {
    if (any(lambda > max(fit$lambda) | lambda < min(fit$lambda))){
      stop('lambda must lie within the range of the fitted coefficient path', call.=FALSE)
    }
    ind <- approx(fit$lambda, seq(fit$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    # linearly interpolate between two lambda
    beta <- (1 - w) * fit$beta[, l, drop = FALSE] + w * fit$beta[, r, drop = FALSE]
    colnames(beta) <- round(lambda, 4)
  } else {  #specify lambda value as index
    beta <- fit$beta[, which, drop = FALSE]
  }
  if (drop == TRUE){
    beta <- drop(beta)
  }
  return(beta)
}

predict.coxkl_highdim <- function(fit, X_new = NULL, lambda = NULL, which = 1:length(fit$lambda),
                                  type = c("link", "response", "vars", "nvars", "groups", "ngroups"), ...){
  beta <- coef.coxkl_highdim(fit, lambda = lambda, which = which, drop = FALSE)
  
  if (type == "vars"){
    return(drop(apply(beta != 0, 2, FUN = which)))
  }
  
  if (type == "groups") {
    if (ncol(beta) == 1) {
      return(unique(fit$group[beta != 0]))
    } else {
      return(drop(apply(beta != 0, 2, function(x) unique(fit$group[x]))))
    }
  }
  
  if (type == "nvars") {
    v <- as.list(apply(beta != 0, 2, FUN = which))
    nvars <- sapply(v, length)
    return(nvars)
  }
  
  if (type == "ngroups") {
    g <- apply(beta != 0, 2, function(x) unique(fit$group[x]))
    ngroups <- sapply(g, length)
    if (length(ngroups) == 1){
      names(ngroups) <- colnames(beta)
    }
    return(ngroups)
  }
  
  # predict response
  if (is.null(data)) {
    stop("Must supply data for predictions", call. = FALSE)
  }
  
  eta <- X_new %*% beta
  if (type == "link") {
    return(eta)
  }
  
  if (type == "response") {
    return(exp(eta))
  }
}



cvf.coxkl_highdim <- function(i, XX, delta, time, stratum, RS, fold, VVH = FALSE, model = c("coxkl_ridge", "coxkl_highdim"), cv.args) {
  cv.args$z <- XX[fold != i, , drop = FALSE]
  cv.args$delta <- delta[fold!=i]
  cv.args$time <- time[fold!=i]
  cv.args$stratum <- stratum[fold!=i]
  cv.args$RS <- RS[fold!=i]
  
  invisible(
    capture.output(
      fit.i <- do.call(model, cv.args)
    )
  )
  nl <- length(fit.i$lambda)
  if (VVH){
    yhat_train <- predict.coxkl_highdim(fit.i, XX[fold != i, , drop = FALSE], type = "link")
    yhat_all <- predict.coxkl_highdim(fit.i, XX, type = "link")
    return(list(nl = nl, yhat_train = yhat_train, yhatall = yhat_all))
  } else {
    yhat_test <- predict.coxkl_highdim(fit.i, XX[fold == i, , drop = FALSE], type = "link")
    return(list(nl = nl, yhat_test = yhat_test))
  }
}


loss.coxkl_highdim <- function(delta, y.hat, stratum, total = TRUE){
  y.hat <- as.matrix(y.hat)
  revCumsum.strat <- function(strat){
    temp.y.hat <- y.hat[which(stratum == strat), , drop=FALSE] 
    temp.rsk <- apply(temp.y.hat, 2, function(x) rev(cumsum(rev(exp(x)))))
    return(temp.rsk)
  }
  rsk <- do.call(rbind, lapply(unique(stratum), function(i) revCumsum.strat(i)))
  
  if (total == TRUE) { #when delta = 0, loss contribution is 0
    return(-2 * (crossprod(delta, y.hat) - crossprod(delta, log(rsk)))) #return a vector (1 * nlambda)
  } else { #when delta = 0, loss contribution is 0
    return(-2 * (y.hat[delta == 1, , drop = FALSE] - log(rsk)[delta == 1, , drop = FALSE])) #return a matrix (n * nlambda)
  }
}



#' Cross-Validation for Lambda in coxkl_highdim (Internal Use)
#'
#' This is an internal helper function used within \code{\link{cv.coxkl_highdim}}. 
#' It performs K-fold cross-validation over the penalty parameter \code{lambda} 
#' for the penalized stratified Cox model, given a fixed \code{eta}.
#'
#' @param z Numeric matrix of predictors (rows = samples, cols = variables).
#' @param delta Event indicator vector (1 = event, 0 = censored).
#' @param time Survival time vector.
#' @param stratum Stratum indicator vector.
#' @param RS External risk score matrix, pre-computed or derived from external \code{beta}.
#' @param eta Numeric, tuning parameter for external information weighting.
#' @param alpha Elastic net mixing parameter (1 = lasso, 0 < alpha < 1 = elastic net).
#' @param nfolds Number of cross-validation folds.
#' @param seed Optional random seed for reproducibility.
#' @param cv.lambda.criteria Character, criteria for lambda selection. 
#'   One of \code{"V&VH"} or \code{"LinPred"}.
#' @param trace.cv Logical; if \code{TRUE}, progress of cross-validation is printed.
#' @param ... Additional arguments passed to \code{\link{coxkl_highdim}}.
#'
#' @return A list containing:
#'   \item{cve}{Vector of cross-validation errors for each lambda.}
#'   \item{beta}{Matrix of fitted coefficients for all lambdas.}
#'   \item{lambda.min}{Selected lambda that minimizes CV error.}
#'   \item{best.beta}{Coefficient estimates corresponding to \code{lambda.min}.}
#'
cv.lambda <- function(
    z, delta, time, stratum, RS, eta, alpha, nfolds, seed = NULL, 
    cv.lambda.criteria = "V&VH", trace.cv = FALSE, ...) {
  cv.lambda.criteria <- match.arg(cv.lambda.criteria, choices = c("V&VH", "LinPred"))
  
  fit.args <- list(...)
  fit.args$z <- z
  fit.args$delta <- delta
  fit.args$time <- time
  fit.args$stratum <- stratum
  fit.args$RS <- RS
  fit.args$beta <- NULL
  fit.args$eta <- eta
  fit.args$alpha <- alpha
  fit.args$returnX <- TRUE
  fit.args$data_sorted <- TRUE
  
  invisible(
    capture.output(
      fit <- do.call("coxkl_highdim", fit.args)
    )
  )
  
  XX <- fit$returnX$XX$std.Z
  delta <- fit$returnX$delta
  time <- fit$returnX$time
  stratum <- fit$returnX$stratum
  RS <- fit$returnX$RS
  
  if (!is.null(seed)) {
    set.seed(seed)
  } else {
    set.seed(NULL)
  }
  fold <- get_fold(nfolds = nfolds, delta = delta, stratum = stratum)
  
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$group <- fit$returnX$XX$g
  cv.args$group.multiplier <- fit$returnX$XX$m
  cv.args$beta <- NULL
  cv.args$alpha <- alpha
  cv.args$eta <- eta
  cv.args$data_sorted <- TRUE
  cv.args$standardize <- FALSE
  
  # Cross-Validation
  Y <- matrix(NA, nrow = length(delta), ncol = length(fit$lambda))
  if (cv.lambda.criteria == "LinPred") {
    for (i in 1:nfolds) {
      if (trace.cv == TRUE){
        cat("Starting CV fold #", i, sep="","...\n")
      }
      res <- cvf.coxkl_highdim(i, XX, delta, time, stratum, RS, fold, VVH = FALSE, model = "coxkl_highdim", cv.args)
      Y[fold == i, 1:res$nl] <- res$yhat_test
    }
    ind <- which(apply(is.finite(Y), 2, all))
    lambda <- fit$lambda[ind]
    L <- loss.coxkl_highdim(delta, Y, stratum, total = FALSE)
    cve <- apply(L, 2, sum)
  } else if (cv.lambda.criteria == "V&VH") {
    cve <- rep(0, length(cv.args$lambda))
    for (i in 1:nfolds) {
      if (trace.cv == TRUE){
        cat("Starting CV fold #", i, sep="","...\n")
      }
      res <- cvf.coxkl_highdim(i, XX, delta, time, stratum, RS, fold, VVH = TRUE, model = "coxkl_highdim", cv.args)
      yhat_train <- res$yhat_train
      yhat_all <- res$yhatall
      L_train <- loss.coxkl_highdim(delta[fold!=i], yhat_train, stratum[fold!=i], total = F)
      L_all <- loss.coxkl_highdim(delta, yhat_all, stratum, total = FALSE)
      cve <- cve + (apply(L_all, 2, sum) - apply(L_train, 2, sum))
    }
  }
  min <- which.min(cve)
  
  best.beta <- fit$beta[, min, drop = FALSE]
  names(cve) <- round(fit$lambda, 4)
  
  return(list(cve = cve,
              beta = fit$beta,
              lambda.min = fit$lambda[min],
              best.beta = best.beta))
}




#' Cross-Validation for Lambda in coxkl_ridge (Internal Use)
#' 
#' Performs K-fold cross-validation over the ridge penalty parameter \code{lambda} 
#' for the penalized stratified Cox model with KL divergence, given a fixed \code{eta}.
#'
#' @param z Numeric matrix of predictors (rows = samples, cols = variables).
#' @param delta Event indicator vector (1 = event, 0 = censored).
#' @param time Survival time vector.
#' @param stratum Stratum indicator vector.
#' @param RS External risk score matrix, pre-computed or derived from external \code{beta}.
#' @param eta Numeric, tuning parameter for external information weighting.
#' @param nfolds Number of cross-validation folds.
#' @param seed Optional random seed for reproducibility.
#' @param cv.lambda.criteria Character, criteria for lambda selection. 
#'   One of \code{"V&VH"} or \code{"LinPred"}.
#' @param trace.cv Logical; if \code{TRUE}, progress of cross-validation is printed.
#' @param ... Additional arguments passed to \code{\link{coxkl_ridge}}.
#'
#' @return A list containing:
#'   \item{cve}{Vector of cross-validation errors for each lambda.}
#'   \item{beta}{Matrix of fitted coefficients for all lambdas.}
#'   \item{lambda.min}{Selected lambda that minimizes CV error.}
#'   \item{best.beta}{Coefficient estimates corresponding to \code{lambda.min}.}
cv.lambda.ridge <- function(
    z, delta, time, stratum, RS, eta, nfolds, seed = NULL, 
    cv.lambda.criteria = "V&VH", trace.cv = FALSE, ...) {
  
  cv.lambda.criteria <- match.arg(cv.lambda.criteria, choices = c("V&VH", "LinPred"))
  
  fit.args <- list(...)
  fit.args$z <- z
  fit.args$delta <- delta
  fit.args$time <- time
  fit.args$stratum <- stratum
  fit.args$RS <- RS
  fit.args$beta <- NULL
  fit.args$eta <- eta
  fit.args$data_sorted <- TRUE
  
  
  invisible(
    capture.output(
      fit <- do.call("coxkl_ridge", fit.args)
    )
  )
  
  if (!is.null(seed)) {
    set.seed(seed)
  } else {
    set.seed(NULL)
  }
  fold <- get_fold(nfolds = nfolds, delta = delta, stratum = stratum)
  
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$beta <- NULL
  cv.args$eta <- eta
  cv.args$data_sorted <- TRUE
  
  # Cross-Validation
  Y <- matrix(NA, nrow = length(delta), ncol = length(fit$lambda))
  if (cv.lambda.criteria == "LinPred") {
    for (i in 1:nfolds) {
      if (trace.cv == TRUE){
        cat("Starting CV fold #", i, sep="","...\n")
      }
      res <- cvf.coxkl_highdim(i, z, delta, time, stratum, RS, fold, VVH = FALSE, model = "coxkl_ridge", cv.args)
      Y[fold == i, 1:res$nl] <- res$yhat_test
    }
    ind <- which(apply(is.finite(Y), 2, all))
    lambda <- fit$lambda[ind]
    L <- loss.coxkl_highdim(delta, Y, stratum, total = FALSE)
    cve <- apply(L, 2, sum)
  } else if (cv.lambda.criteria == "V&VH") {
    cve <- rep(0, length(cv.args$lambda))
    for (i in 1:nfolds) {
      if (trace.cv == TRUE){
        cat("Starting CV fold #", i, sep="","...\n")
      }
      res <- cvf.coxkl_highdim(i, z, delta, time, stratum, RS, fold, VVH = TRUE, model = "coxkl_ridge", cv.args)
      yhat_train <- res$yhat_train
      yhat_all <- res$yhatall
      L_train <- loss.coxkl_highdim(delta[fold!=i], yhat_train, stratum[fold!=i], total = F)
      L_all <- loss.coxkl_highdim(delta, yhat_all, stratum, total = FALSE)
      cve <- cve + (apply(L_all, 2, sum) - apply(L_train, 2, sum))
    }
  }
  min <- which.min(cve)
  
  best.beta <- fit$beta[, min, drop = FALSE]
  names(cve) <- round(fit$lambda, 4)
  
  return(list(cve = cve,
              beta = fit$beta,
              lambda.min = fit$lambda[min],
              best.beta = best.beta))
}

