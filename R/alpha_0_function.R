#' alpha_0_function
#'
#' @param Y_T Data
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return A random starting value of alpha is generated
#' @export
#'
#' @examples alpha_0_function(Y_T,N,K)
alpha_0_function = function(Y_T,N,K){
  alpha = matrix(0,nrow = N,ncol = K+2)
  for(n in 1:N){
    alpha[n,1] = runif(1,min(Y_T),max(Y_T))
    alpha[n,2:(K+1)] = rnorm(K,0,0.5)
  }
  alpha[,K+2] = runif(N,0,sd(Y_T)/2) #Ist dennoch (sigma_u)^2
  return(alpha)
}
