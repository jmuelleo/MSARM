
#' Example4_alpha_next_function_lastbeta
#' Example4_alpha_next_function_lastbeta allows its user to estimate the alpha matrix of the next iteration for Example 4 setups, based on Müller (2025, page 19-21)
#' @param Y_T Data
#' @param alpha alpha matrix from the last iteration
#' @param P P matrix from the last iteration
#' @param Switcher Switching vector which indicates which parameter switch
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Number of time periods the conditional log-likelihood conditions on
#'
#' @return Returns a matrix of the dimensions (N x K+2)
#' @export
#'
#' @examples Example4_alpha_next_function_lastbeta(Y_T,alpha,P,Switcher,N,K,m)
Example4_alpha_next_function_lastbeta = function(Y_T,alpha,P,Switcher,N,K,m){
  T = length(Y_T)
  zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K) #smoothed inference over the regimes
  p = zeta_YT[,(K+K+1+m):T]

  X = matrix(0,ncol = K+1,nrow = T-K-K-m)

  Y = Y_T[(K+K+1+m):T]

  X[,1] = rep(1,T-K-K-m)
  for(k in 1:K){
    X[,k+1] = Y_T[(K+K+1+m-k):(T-k)] #Lags of the time series are created
  }

  coef_Switcher = Switcher[-(K+2)]

  X_F = if(sum(!(Switcher[-(K+2)])) > 0){
    as.matrix(X[,!coef_Switcher])}else{ #Creates the data matrix for the coefficients that are fixed
      c(0)}
  X_S = as.matrix(X[, coef_Switcher]) #creates the data matrix for the coefficients that switch
  F_count = sum(coef_Switcher == FALSE) #number of fixed parameters
  S_count = sum(coef_Switcher == TRUE) #number of switching parameters

  Beta_F_value = Example4_Beta_F_function(Y_T,alpha,P,Switcher,N,K,m) #creates the fixed beta based on the data from the previous iteration, based on the formulia in Müller (2025, page 21)
  Beta = matrix(5,ncol = K+1, nrow = N)
  Residuals = matrix(5,ncol = N,nrow = T-K-K-m)

  #creates the total beta vector, including the beta_F and beta_S vector
  for(i in 1:N){

    Y_Star = Y*sqrt(p[i,])
    X_S_Star = X_S*c(sqrt(p[i,]))
    X_F_Star = X_F*c(sqrt(p[i,]))

    Y_Star_adj = Y_Star - X_F_Star%*%as.matrix(Beta_F_value)

    erg = lm(Y_Star_adj ~ X_S_Star-1) #creates the beta_S vector, based on Müller (2025, page 21)

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

  Sigma2 = rep(sum(apply((Residuals^2)/(T-K-K-m),2,sum)),N) #estimates the error term variance as described in Müller (2025, page 21)

  alpha = cbind(Beta,Sigma2)
  return(alpha)
}

