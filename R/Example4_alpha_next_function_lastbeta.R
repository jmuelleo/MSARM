#' Title
#'
#' @param Y_T Data
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param P Transpose of the transition matrix Pi
#' @param Switcher Switching-vector, a vector of (TRUE/FALSE) of the length K+2, which indicates which parameters are supposed to switch
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return Computes, for Example 4 setups (Subset of the coefficients switches, error term variance is fixed) the alpha of the next iteration, based on the results from the last iteration
#' @export
#'
#' @examples Example4_alpha_next_function_lastbeta(Y_T,alpha,P,Switcher,N,K)
Example4_alpha_next_function_lastbeta = function(Y_T,alpha,P,Switcher,N,K,m){
  T = length(Y_T)
  zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K)
  p = zeta_YT[,(K+K+1+m):T]

  X = matrix(0,ncol = K+1,nrow = T-K-K-m)

  Y = Y_T[(K+K+1+m):T]

  X[,1] = rep(1,T-K-K-m)
  for(k in 1:K){
    X[,k+1] = Y_T[(K+K+1+m-k):(T-k)]
  }

  coef_Switcher = Switcher[-(K+2)]

  X_F = if(sum(!(Switcher[-(K+2)])) > 0){
    as.matrix(X[,!coef_Switcher])}else{
      c(0)}
  X_S = as.matrix(X[, coef_Switcher])
  F_count = sum(coef_Switcher == FALSE)
  S_count = sum(coef_Switcher == TRUE)

  Beta_F_value = Example4_Beta_F_function(Y_T,alpha,P,Switcher,N,K,m)
  Beta = matrix(5,ncol = K+1, nrow = N)
  Residuals = matrix(5,ncol = N,nrow = T-K-K-m)

  for(i in 1:N){

    Y_Star = Y*sqrt(p[i,])
    X_S_Star = X_S*c(sqrt(p[i,]))
    X_F_Star = X_F*c(sqrt(p[i,]))

    Y_Star_adj = Y_Star - X_F_Star%*%as.matrix(Beta_F_value)

    erg = lm(Y_Star_adj ~ X_S_Star-1)

    Beta_S  = coef(erg)

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



    Residuals[,i] = erg$residuals
  }

  Sigma2 = rep(sum(apply((Residuals^2)/(T-K-K-m),2,sum)),N)

  alpha = cbind(Beta,Sigma2)
  return(alpha)
}

