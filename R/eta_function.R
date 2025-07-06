#' eta_function
#' The eta-vector, defined as in Hamilton (1994, page 691) is computed
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return returns a matrix of the dimensions (T x N) containing the values of the conditional density function
#' @export
#'
#' @examples eta_function(alpha,Y_T,N,K)
eta_function = function(alpha,Y_T,N,K){
  T = length(Y_T)
  eta = matrix(5,ncol = T,nrow = N)
  for(t in (K+1):T){ #presample data is necessary, therefore + K -> (total time series has the length T = l + K, core time series has the length l)
    for(h in 1:N){
      eta[h,t] = dnorm(Y_T[t], #follows the eta-vector has defined in Hamilton (1994, page 691)
                       mean = c(1,Y_T[(t-1):(t-K)])%*%alpha[h,-(K+2)],
                       sd = sqrt(alpha[h,(K+2)]))

    }
  }
  return(eta + 0.000000001) #At the beginning of the algorithm extreme values are possible, therefore a very small number is added to avoid computation error
}
