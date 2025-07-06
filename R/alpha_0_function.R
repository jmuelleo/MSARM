
#' alpha_0_function
#' A random starting value of alpha is generated so that the EM algorithm can start
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return
#' @export
#'
#' @examples alpha_0_function = function(Y_T,N,K)
#'
alpha_0_function = function(Y_T,N,K){
  T = length(Y_T)
  alpha = matrix(0,nrow = N,ncol = K+2)
  for(n in 1:N){
    alpha[n,1] = runif(1,min(Y_T),max(Y_T)) #intercept is drawn from a uniform distribution ranging from the minimum to the maximum of the time series
    alpha[n,2:(K+1)] = rnorm(K,0,0.5) #the coefficients are drawn from a normal distribution with mean 0 and standard deviation 0.5
  }
  alpha[,K+2] = runif(N,0,sd(Y_T)/2) #the error term variance is first set between 0 and the standard deviation of the total time series divided by 2
  return(alpha)
}
