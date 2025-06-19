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
#' @examples Example5_alpha_next_function_lastbeta(Y_T,alpha,P,Switcher,N,K)
Example5_alpha_next_function_lastbeta = function(Y_T,alpha,P,Switcher,N,K){

  zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K)
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

  Beta_F_value = Example5_Beta_F_function(Y_T,alpha,P,Switcher,N,K)
  sigma = Example5_sigma_function(Y_T,alpha,P,Switcher,N,K)
  Beta = matrix(5,ncol = K+1, nrow = N)
  Residuals = matrix(5,ncol = N,nrow = T-K-K-m)

  for(i in 1:N){

    Y_Star = Y*(sqrt(p[i,])/sigma[i])
    X_S_Star = X_S*(sqrt(p[i,])/sigma[i])
    X_F_Star = X_F*(sqrt(p[i,])/sigma[i])

    Y_Star_adj = Y_Star - X_F_Star%*%Beta_F_value

    if(dim(X_S_Star)[2] > 0){
      erg = lm(Y_Star_adj ~ X_S_Star-1)
    }else{erg = 0}

    #erg = lm(Y_Star_adj ~ X_S_Star-1)

    if(dim(X_S_Star)[2] > 0){
      Beta_S  = coef(erg)}else{
        Beta_S  = 0
      }

    result = rep(0,K+1)
    S_i = 1
    F_i = 1
    for(k in 1:(K+1)){
      if(Switcher[k] == TRUE){
        result[k] = Beta_S[S_i]
        S_i = S_i + 1
      }else{
        result[k] = Beta_F_value[F_i]
        F_i = F_i + 1
      }
    }
    Beta[i,] = result
  }

  alpha = cbind(Beta,sigma^2)
  return(alpha)
}

