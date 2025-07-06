
#' Example5_Beta_F_function
#' Example5_Beta_F_function allows its user to compute beta_F as part of the EM algorithm for Example 5 setups, for details see Müller (2025, page 23-24)
#' @param Y_T Data
#' @param alpha alpha matrix from the previous iteration
#' @param P P matrix from the previous iteration
#' @param Switcher Switching vector indicating whether a parameter is fixed or switches
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Number of time periods the conditional likelihood conditions on
#'
#' @return Returns a vector representing the estimates for beta_F based on the results from the last iteration
#' @export
#'
#' @examples Example5_Beta_F_function(Y_T,alpha,P,Switcher,N,K,m)
Example5_Beta_F_function = function(Y_T,alpha,P,Switcher,N,K,m){
  T = length(Y_T)
  zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K) #smoothed inference regarding the regimes
  p = zeta_YT[,(K+K+1):T]

  X = matrix(0,ncol = K+1,nrow = T-K-K)

  Y = Y_T[(K+K+1):T]

  X[,1] = rep(1,T-K-K)
  for(k in 1:K){
    X[,k+1] = Y_T[(K+K+1-k):(T-k)] #Creates the Lags of the time series
  }

  coef_Switcher = Switcher[-(K+2)]

  X_F = as.matrix(X[,!coef_Switcher]) #Creates the data matrix for the fixed coefficients
  X_S = as.matrix(X[, coef_Switcher]) #Creates the data matrix for the switching coefficients

  F_count = sum(coef_Switcher == FALSE) #Number of fixed parameters
  S_count = sum(coef_Switcher == TRUE) #Number of switching parameters


  sigma = Example5_sigma_function(Y_T,alpha,P,Switcher,N,K,m) #computes the error term standard deviation based on the results from the previous iteration and the formula presented in Müller (2025, page 24)


  #creates the divisor matrix as the first part of beta_F, based on the fomula in Müller (2025, page 24)
  start_matrix = matrix(0,nrow = F_count,ncol = F_count)
  for(t in (m+1):(T-K-K)){
    xF_t = X_F[t,]
    matrix_t = xF_t%*%t(xF_t)

    multiplier_vector = rep(0,N)
    for(j in 1:N){
      multiplier_vector[j] = p[j,t]/((sigma[j])^2)
    }
    multiplier = sum(multiplier_vector)

    start_matrix = start_matrix + (matrix_t*multiplier)
  }
  divisor_matrix = start_matrix


  #creates the second part of the equation for beta_F, based on the formula in Müller (2025, page 24)
  start_vector = matrix(0,ncol = 1,nrow = F_count)
  for(t in (m+1):(T-K-K)){
    for(j in 1:N){
      y_t = Y[t]
      xF_t = X_F[t,]
      xS_t = X_S[t,]
      betaS_j = (alpha[j,c(coef_Switcher,FALSE)])

      vector_t = ((y_t - (xS_t%*%betaS_j))%*%xF_t)*(p[j,t]/(sigma[j]^2))
      start_vector = start_vector + t(vector_t)
    }
  }
  upper_vector = start_vector

  #Computes the total beta_F vector
  Beta_F = solve(divisor_matrix)%*%upper_vector
  return(Beta_F)
}
