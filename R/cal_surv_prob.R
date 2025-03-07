cal_surv_prob <- function(z, delta, time, beta){
  delta = delta[order(time)]
  z = z[order(time),]
  time = time[order(time)]
  
  z_mat <- as.matrix(z)
  delta_mat <- as.matrix(delta)
  beta <- as.matrix(beta)
  diff = ddloglik_S0(z_mat,delta_mat,beta)
  S0   = diff$S0
  
  Lambda0=cumsum(delta/S0)
  tmax = length(Lambda0)
  n = nrow(z_mat)
  S <- matrix(rep(0, (n*tmax)), n, tmax)
  for (i in 1:tmax){
    S[,i]=exp(-Lambda0[i]*exp(z_mat%*%beta))
  }
  
  return_list <- list("S"=S)
  return(return_list)
}


#' Calculate the likelihood
#' 
#' @export
loss_fn <- function(z, delta, time, beta){
  delta = delta[order(time)]
  z = z[order(time),]
  time = time[order(time)]
  
  z_mat <- as.matrix(z)
  delta_mat <- as.matrix(delta)
  beta <- as.matrix(beta)
  diff = loss_fn_cpp(z_mat,delta_mat,beta)
  S0   = diff$loglik

  return(S0)
}

