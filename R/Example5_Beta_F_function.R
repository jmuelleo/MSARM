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
#' @examples Example5_Beta_F_function(Y_T,alpha,P,Switcher,N,K)
Example5_Beta_F_function = function(Y_T,alpha,P,Switcher,N,K){

  zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K)
  p = zeta_YT[,(K+K+1):T]

  X = matrix(0,ncol = K+1,nrow = T-K-K)

  Y = Y_T[(K+K+1):T]

  X[,1] = rep(1,T-K-K)
  for(k in 1:K){
    X[,k+1] = Y_T[(K+K+1-k):(T-k)]
  }

  coef_Switcher = Switcher[-(K+2)]

  X_F = as.matrix(X[,!coef_Switcher])
  X_S = as.matrix(X[, coef_Switcher])

  F_count = sum(coef_Switcher == FALSE)
  S_count = sum(coef_Switcher == TRUE)


  sigma = Example5_sigma_function(Y_T,alpha,P,Switcher,N,K)


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

  Beta_F = solve(divisor_matrix)%*%upper_vector
  return(Beta_F)
}
