#' Title
#'
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param P Transpose of the transition matrix Pi
#' @param Y_T Data
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return Evaluates the conditional log-likelihood function
#' @export
#'
#' @examples Log_Likelihood_function(alpha,P,Y_T,N,K)
Log_Likelihood_function = function(alpha,P,Y_T,N,K){
  zeta_list = zeta_Yt_function(alpha,P,Y_T,N,K)
  zeta_t_t = zeta_list[[1]]
  zeta_t1_t = zeta_list[[2]]
  eta = eta_function(alpha,Y_T,N,K)

  f_sammler = rep(5,T-K)
  for(t in (K+1):T){
    f_sammler[t-K] = sum(zeta_t1_t[,t-1]*eta[,t])
  }

  Log_Likelihood_value = sum(log(f_sammler + 0.000000001)) #damit durch Rundungen kein echten Nullen enthalten sind
  return(Log_Likelihood_value)
}
