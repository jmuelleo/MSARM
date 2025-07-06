
#' zeta_Yt_function
#' zeta_Yt_function allows its user to get optimal inference over the regimes as described in Hamilton (1994, page 692)
#' @param alpha alpha matrix from the last iteration
#' @param P P matrix from the last iteration
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return Returns the optimal inference over the regimes
#' @export
#'
#' @examples zeta_Yt_function(alpha,P,Y_T,N,K)
zeta_Yt_function = function(alpha,P,Y_T,N,K){
  T = length(Y_T)
  zeta_1_0 = zeta_1_0_function(Y_T,N,K) #random starting point for the optimal inference over the regimes
  eta = eta_function(alpha = alpha ,Y_T = Y_T,N = N,K = K) #eta vector as described in Hamilton (1994, page 691)

  zeta_t_t = matrix(5,ncol = T,nrow = N)
  zeta_t1_t = matrix(5,ncol = T,nrow = N)

  zeta_t1_t[,K] = zeta_1_0

  zeta_t_t[,K+1] = (zeta_t1_t[,K]*eta[,K+1])/sum(zeta_t1_t[,K]*eta[,K+1])
  zeta_t1_t[,K+1] = P%*%zeta_t_t[,K+1]
  #zeta_t1_t contains the t+1 value in column t! "Given/Conditional on" determines the column, t given -1 is column t-1

  for(t in (K+2):T){
    zeta_t_t[,t] = (zeta_t1_t[,t-1]*eta[,t])/sum(zeta_t1_t[,t-1]*eta[,t]) #computation based on the formulas in Hamilton (1994, page 692)
    zeta_t1_t[,t] = P%*%zeta_t_t[,t]
  }

  return(list(tt = zeta_t_t, t1t =  zeta_t1_t))

}
