# 
# loss.MIGroupLasso <- function(y, eta, K, total=TRUE) {
#   l <- dim(eta)[2]
#   ind <- order(y[,1])
#   d <- as.integer(y[ind,2])
#   N <- length(d)
#   ind2 <- NULL
#   for(i in 1:K){
#     ind2 <- c(ind2, ind + N*(i-1))
#   }
#   
#   eta <- eta[ind2, , drop=FALSE]
#   r   <- matrix(0,N*K,l)
#   for(i in 1:K){
#     samplei <- (N*(i-1)+1) : (N*i)
#     r[samplei,] <- apply(eta[samplei,], 2, function(x) rev(cumsum(rev(exp(x)))))
#   }
#   
#   loglik = 0
#   for (i in 1:K) {
#     samplei  <- (N*(i-1)+1) : (N*i) & d==1
#     loglik = loglik -2*(eta[samplei, , drop=FALSE] - log(r)[samplei, , drop=FALSE])
#   }
#   
#   return(loglik)
# }

loss.MIGroupLasso <- function(y, eta, K, total=TRUE) {
  l <- dim(eta)[2]
  ind <- order(y[,1])
  #print(ind)
  d <- as.integer(y[ind,2])
  N <- length(d)
  ind2 <- NULL
  for(i in 1:K){
    ind2 <- c(ind2, ind + N*(i-1))
  }
  
  eta <- eta[ind2, , drop=FALSE]
  r   <- matrix(0,N*K,l)
  loglik = 0
  for(i in 1:K){
    samplei <- (N*(i-1)+1) : (N*i)
    #r[samplei,] <- apply(eta[samplei,], 2, function(x) rev(cumsum(rev(exp(x)))))
    r  <- apply(eta[samplei,], 2, function(x) rev(cumsum(rev(exp(x)))))
    indexi_death  <- which((N*(i-1)+1) : (N*i) & d==1)
    #loglik = loglik - 2*(eta[indexi_death, , drop=FALSE] - log(r)[d==1, , drop=FALSE])
    
    
    indexi_death  <- (N*(i-1)+1) : (N*i)
    loglik = loglik - 2*(crossprod(d, eta[indexi_death, , drop=FALSE]) - crossprod(d, log(r)))
  }
  
  
  
  return(loglik)
}




