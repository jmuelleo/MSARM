
#' @title Example3_alpha_next_function
#' @description The next alpha for Example 3 setups (all parameters switch) is computed as part of the EM algorithm
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param P transpose of the transition matrix (of the Markov-Chain) from the last iteration
#' @param Y_T Data
#' @param N Number of regimes
#' @param K Lag-order of the AR model
#' @param m Number of time periods the conditional likelihood conditions on
#'
#' @return returns a matrix of the dimensions (N x K+2)
#' @export
#'
#' @examples Example3_alpha_next_function(alpha,P,Y_T,N,K,m)
Example3_alpha_next_function = function(alpha,P,Y_T,N,K,m){
  T = length(Y_T)
  zeta_list = zeta_Yt_function(alpha,P,Y_T,N,K) #optimal inference regarding the regimes
  zeta_t_t = zeta_list[[1]]
  zeta_t1_t = zeta_list[[2]]
  zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K) #smoothed inference regarding the regimes
  zeta = zeta_YT[,(K+1+K):T] #presample data is omitted

  X = matrix(0,ncol = K,nrow = T-K-K)
  Beta = matrix(0,ncol = K+1,nrow = N)
  Sigma = matrix(5,ncol = 1,nrow = N)

  Y = Y_T[(K+1+K):T] #presample data is omiitted

  for(k in 1:(K)){
    X[,k] = (Y_T[(K+1+K-k):(T-k)]) #Lags of the time series are created
  }

  sigma_j = function(beta_j,j,Y_T,N,K){
    sqrt(sum((Y-cbind(rep(1,T-K-K),X)%*%beta_j)^2*zeta[j,])/sum(zeta[j,])) #sigma is computed based on Müller (2025, page 19), with the approximation that
  # the values from the last iteration are utilized to compute sigma_j
  }

  for(i in 1:N){
    sigma = sigma_j(beta_j = alpha[i,-(K+2)],j = i,Y_T = Y_T,N = N,K = K)
    result = lm(Y ~ X, weights = (zeta_YT[i,(K+1+K):T])/(sigma^2)) #weighted regression is conducted based on Müller (2025, page 18-19)
    Beta[i,] = result$coefficients
    Sigma[i,] = (sigma_j(Beta[i,],i,Y_T,N,K))^2
  }



  alpha = cbind(Beta,Sigma)
  return(alpha)
}
