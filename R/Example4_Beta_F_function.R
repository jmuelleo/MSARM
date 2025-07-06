
#' Example4_Beta_F_function
#' Example4_Beta_F_function allows its user to estimate the fixed coefficients vector called beta_F, based on the results from the previous iteration and the formulas in Müller (2025, page 21)
#' @param Y_T Data
#' @param alpha alpha matrix from the previous iteration
#' @param P P matrix from the previous iteration
#' @param Switcher Switching vector indicating which parameters switch
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Number of time periods the conditional log-likelihood conditions on
#'
#' @return Returns a vector representing beta_F based on the results from the last iteration
#' @export
#'
#' @examples Example4_Beta_F_function(Y_T,alpha,P,Switcher = Switcher,N,K,m)
Example4_Beta_F_function = function(Y_T,alpha,P,Switcher = Switcher,N,K,m){
  if(sum(!(Switcher[-(K+2)])) > 0){
    zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K) #smoothed inference over the regimes
    T = length(Y_T)
    X = matrix(0,ncol = K+1,nrow = T-K-K)

    Y = Y_T[(K+K+1):T]

    X[,1] = rep(1,T-K-K)
    for(k in 1:K){
      X[,k+1] = Y_T[(K+K+1-k):(T-k)] #Creating the Lags of the time series
    }

    coef_Switcher = Switcher[-(K+2)]

    X_F = as.matrix(X[,!coef_Switcher])
    X_S = as.matrix(X[, coef_Switcher])

    F_count = sum(coef_Switcher == FALSE) #Number of fixed parameters
    S_count = sum(coef_Switcher == TRUE) #Number of switching parameters


    #create the divisor matrix as in Müller (2025, page 21) for beta_F
    start_matrix = matrix(0,nrow = F_count,ncol = F_count)
    for(t in (m+1):(T-K-K)){
      xF_t = X_F[t,]
      matrix_t = xF_t%*%t(xF_t)
      start_matrix = start_matrix + matrix_t
    }
    divisor_matrix = start_matrix


    #creates the second part of the formula for beta_F
    start_vector = matrix(0,ncol = 1,nrow = F_count)
    for(t in (m+1):(T-K-K)){
      for(j in 1:N){
        y_t = Y[t]
        xF_t = X_F[t,]
        xS_t = X_S[t,]
        betaS_j = (alpha[j,c(coef_Switcher,FALSE)])
        p = zeta_YT[,(K+K+1):T]

        vector_t = ((y_t - (xS_t%*%betaS_j))%*%xF_t)*p[j,t]
        start_vector = start_vector + t(vector_t)
      }
    }
    upper_vector = start_vector

    Beta_F = solve(divisor_matrix)%*%upper_vector #compute the total beta_F vector
    return(Beta_F)
  }else{
    return(c(0))
  }
}
