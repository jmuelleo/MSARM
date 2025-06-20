#' Title
#'
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param P Transpose of the transition matrix Pi
#' @param Y_T Data
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return Creates the smoothed inference for the regimes,utilizing the Kim-Algorithm.
#'  The theory can be found in Kim(1993) or Hamilton(1994) p.694 and p.699-700
#' @export
#'
#' @examples zeta_YT_function(alpha,P,Y_T,N,K)
zeta_YT_function = function(alpha,P,Y_T,N,K){
  T = length(Y_T)
  zeta_list = zeta_Yt_function(alpha,P,Y_T,N,K)
  zeta_t_t = zeta_list[[1]]
  zeta_t1_t = zeta_list[[2]]

  zeta_t_T = matrix(5,nrow = N,ncol = T)
  zeta_t_T[,T] = zeta_t_t[,T]

  for(t in (T-1):(K+1)){
    zeta_t_T[,t] = zeta_t_t[,t] * (t(P)%*%(zeta_t_T[,t+1]/zeta_t1_t[,t]))
  }
  return(zeta_t_T)
}
