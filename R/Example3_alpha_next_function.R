#' Title
#'
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param P Transpose of the transition matrix Pi
#' @param Y_T Data
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return Computes,for Example 3 setups (All Switching), the alpha matrix for the next iteration, based on the last iteration.
#' @export
#'
#' @examples Example3_alpha_next_function(alpha,P,Y_T,N,K)
Example3_alpha_next_function = function(alpha,P,Y_T,N,K,m){
  T = length(Y_T)
  zeta_list = zeta_Yt_function(alpha,P,Y_T,N,K)
  zeta_t_t = zeta_list[[1]]
  zeta_t1_t = zeta_list[[2]]
  zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K)
  zeta = zeta_YT[,(K+1+K):T]

  X = matrix(0,ncol = K,nrow = T-K-K)
  Beta = matrix(0,ncol = K+1,nrow = N)
  Sigma = matrix(5,ncol = 1,nrow = N)

  Y = Y_T[(K+1+K):T]

  for(k in 1:(K)){
    X[,k] = (Y_T[(K+1+K-k):(T-k)])
  }

  sigma_j = function(beta_j,j,Y_T,N,K){
    sqrt(sum((Y-cbind(rep(1,T-K-K),X)%*%beta_j)^2*zeta[j,])/sum(zeta[j,]))
  }

  for(i in 1:N){
    sigma = sigma_j(beta_j = alpha[i,-(K+2)],j = i,Y_T = Y_T,N = N,K = K)
    result = lm(Y ~ X, weights = (zeta_YT[i,(K+1+K):T])/(sigma^2))
    Beta[i,] = result$coefficients
    Sigma[i,] = (sigma_j(Beta[i,],i,Y_T,N,K))^2
  }



  alpha = cbind(Beta,Sigma)
  return(alpha)
}
