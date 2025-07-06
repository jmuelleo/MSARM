
#' zeta_YT_function
#' zeta_YT_function allows its user to get smoothed inference over the regimes as described in Hamilton (1994, page 694)
#' @param alpha alpha matrix from the last iteration
#' @param P P matrix from the last iteration
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return Returns the smoothed inference over the regimes
#' @export
#'
#' @examples zeta_YT_function(alpha,P,Y_T,N,K)
zeta_YT_function = function(alpha,P,Y_T,N,K){
  T = length(Y_T)
  zeta_list = zeta_Yt_function(alpha,P,Y_T,N,K) #optimal inference over the regimes
  zeta_t_t = zeta_list[[1]]
  zeta_t1_t = zeta_list[[2]]

  zeta_t_T = matrix(5,nrow = N,ncol = T)
  zeta_t_T[,T] = zeta_t_t[,T]

  for(t in (T-1):(K+1)){
    zeta_t_T[,t] = zeta_t_t[,t] * (t(P)%*%(zeta_t_T[,t+1]/zeta_t1_t[,t])) #Estimation based on the formulas shown in Hamilton (1994, page 694)
  }
  return(zeta_t_T)
}
