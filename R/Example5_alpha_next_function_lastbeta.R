
#' Example5_alpha_next_function_lastbeta
#' Example5_alpha_next_function_lastbeta allows its user to estimate the alpha matrix of the next iteration as part of the EM algorithm for setups of the form of Example 5. A more detailed explanation can be found in Müller (2025, page 21-24)
#' @param Y_T Data
#' @param alpha alpha matrix from the previous iteration
#' @param P P matrix from the previous iteration
#' @param Switcher Switching vector indicating which parameters switch
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Number of time periods the conditional likelihood conditions on
#'
#' @return Returns a matrix of the dimensions (N x K+2)
#' @export
#'
#' @examples Example5_alpha_next_function_lastbeta(Y_T,alpha,P,Switcher,N,K,m)
Example5_alpha_next_function_lastbeta = function(Y_T,alpha,P,Switcher,N,K,m){
  T = length(Y_T)
  zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K) #smoothed inference over the regimes
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

  Beta_F_value = Example5_Beta_F_function(Y_T,alpha,P,Switcher,N,K,m) #Estimate the beta_F value based on the results from the previous iteration and the formula in Müller (2025, page 23-24)
  sigma = Example5_sigma_function(Y_T,alpha,P,Switcher,N,K,m) #Estimate the error term standard deviation based on the results from the previous iteration and the formula in Müller (2025, page 23-24)
  Beta = matrix(5,ncol = K+1, nrow = N)
  Residuals = matrix(5,ncol = N,nrow = T-K-K-m)

  #value of the beta vector is created, inclduing both beta_F and beta_S.
  for(i in 1:N){

    Y_Star = Y*(sqrt(p[i,])/sigma[i])
    X_S_Star = X_S*(sqrt(p[i,])/sigma[i])
    X_F_Star = X_F*(sqrt(p[i,])/sigma[i])

    Y_Star_adj = Y_Star - X_F_Star%*%Beta_F_value

    if(dim(X_S_Star)[2] > 0){
      erg = lm(Y_Star_adj ~ X_S_Star-1) #Creating the beta_S vector, based on Müller (2025, page 22)
    }else{erg = 0}

    #erg = lm(Y_Star_adj ~ X_S_Star-1)

    if(dim(X_S_Star)[2] > 0){
      Beta_S  = coef(erg)}else{
        Beta_S  = 0
      }

    #bringing the beta_F and beta_S vector together in the beta vector
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

