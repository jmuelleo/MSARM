#' Title
#'
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param P Transpose of the transition matrix Pi
#' @param Y_T Data
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return Estimates the optimal inference (Bayes Formula) about regimes
#' @export
#'
#' @examples zeta_Yt_function(alpha,P,Y_T,N,K)
zeta_Yt_function = function(alpha,P,Y_T,N,K){
  T = length(Y_T)
  zeta_1_0 = zeta_1_0_function(Y_T,N,K)
  eta = eta_function(alpha = alpha ,Y_T = Y_T,N = N,K = K)

  zeta_t_t = matrix(5,ncol = T,nrow = N)
  zeta_t1_t = matrix(5,ncol = T,nrow = N)

  zeta_t1_t[,K] = zeta_1_0

  zeta_t_t[,K+1] = (zeta_t1_t[,K]*eta[,K+1])/sum(zeta_t1_t[,K]*eta[,K+1])
  zeta_t1_t[,K+1] = P%*%zeta_t_t[,K+1]
  #zeta_t1_t beinhaltet den t+1 Wert in der t Spalte sozusagen! "GEGEBEN" bestimmt die spalte, t gegeben -1 ist in spalte t-1

  for(t in (K+2):T){
    zeta_t_t[,t] = (zeta_t1_t[,t-1]*eta[,t])/sum(zeta_t1_t[,t-1]*eta[,t])
    zeta_t1_t[,t] = P%*%zeta_t_t[,t]
  }

  return(list(tt = zeta_t_t, t1t =  zeta_t1_t))

}
