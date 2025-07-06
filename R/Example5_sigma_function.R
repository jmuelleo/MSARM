
#' Example5_sigma_function
#' Example5_sigma_function allows its user to esimate the error term standard deviation for Example 5 setups, based on the results from previous iterations and the formula presented in Müller (2025, page 24)
#' @param Y_T Data
#' @param alpha alpha matrix from the previous iteration
#' @param P P matrix from the previous iteration
#' @param Switcher Switching vector indicating whether a parameter switches or not
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Number of time periods the conditional likelihood conditions on
#'
#' @return Returns a vector representing the estimates for the standard deviation of the error term, based on the results from the last iteration
#' @export
#'
#' @examples Example5_sigma_function(Y_T,alpha,P,Switcher,N,K,m)
Example5_sigma_function = function(Y_T,alpha,P,Switcher,N,K,m){

  T = length(Y_T)
  zeta_YT = zeta_YT_function(alpha = alpha,P = P ,Y_T = Y_T,N = N,K = K) #smoothed inference for the regimes
  p = zeta_YT[,(K+K+1+m):T]

  X = matrix(0,ncol = K+1,nrow = T-K-K-m)

  Y = Y_T[(K+K+1+m):T]

  X[,1] = rep(1,T-K-K-m)
  for(k in 1:K){
    X[,k+1] = Y_T[(K+K+1+m-k):(T-k)] #Creates the Lags of the time series
  }

  coef_Switcher = Switcher[-(K+2)]

  X_F = as.matrix(X[,!coef_Switcher]) #Creates the data matrix for the fixed coefficients
  X_S = as.matrix(X[, coef_Switcher]) #Creates the data matrix for the switching coefficients
  F_count = sum(coef_Switcher == FALSE) #Number of fixed parameters
  S_count = sum(coef_Switcher == TRUE) #Number of switching parameters

  sigma = rep(0,N)
  for(j in 1:N){

    upper_j = sum(((Y - X%*%alpha[j,-(K+2)])^2)*p[j,]) #creates the numerator of the formula presented in Müller (2025, page 24)
    lower_j = sum(p[j,]) #creates the denominator of the formula presented in Müller (2025, page 24)

    sigma[j] = sqrt(upper_j/lower_j) #compute sigma_j

  }

  return(sigma)
}
