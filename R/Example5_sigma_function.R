#' Title
#'
#' @param Y_T Data
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param P Transpose of the transition matrix Pi
#' @param Switcher Switching-vector, a vector of (TRUE/FALSE) of the length K+2, which indicates which parameters are supposed to switch
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return
#' @export
#'
#' @examples Example5_sigma_function(Y_T,alpha,P,Switcher,N,K)
Example5_sigma_function = function(Y_T,alpha,P,Switcher,N,K,m){

  T = length(Y_T)
  zeta_YT = zeta_YT_function(alpha = alpha,P = P ,Y_T = Y_T,N = N,K = K)
  p = zeta_YT[,(K+K+1+m):T]

  X = matrix(0,ncol = K+1,nrow = T-K-K-m)

  Y = Y_T[(K+K+1+m):T]

  X[,1] = rep(1,T-K-K-m)
  for(k in 1:K){
    X[,k+1] = Y_T[(K+K+1+m-k):(T-k)]
  }

  coef_Switcher = Switcher[-(K+2)]

  X_F = as.matrix(X[,!coef_Switcher])
  X_S = as.matrix(X[, coef_Switcher])
  F_count = sum(coef_Switcher == FALSE)
  S_count = sum(coef_Switcher == TRUE)

  sigma = rep(0,N)
  for(j in 1:N){

    upper_j = sum(((Y - X%*%alpha[j,-(K+2)])^2)*p[j,])
    lower_j = sum(p[j,])

    sigma[j] = sqrt(upper_j/lower_j)

  }

  return(sigma)
}
