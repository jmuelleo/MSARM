
#' Log_Likelihood_function
#' Log_Likelihood_function allows its user to compute the value of the log-likelihood function based on Hamilton (1994, page 692)
#' @param alpha alpha vector from the previous iteration
#' @param P P matrix from the previous iteration
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return
#' @export
#'
#' @examples Log_Likelihood_function(alpha,P,Y_T,N,K)
Log_Likelihood_function = function(alpha,P,Y_T,N,K){
  T = length(Y_T)
  zeta_list = zeta_Yt_function(alpha,P,Y_T,N,K) #optimal inference over the regimes
  zeta_t_t = zeta_list[[1]]
  zeta_t1_t = zeta_list[[2]]
  eta = eta_function(alpha,Y_T,N,K) #eta vector as described in Hamilton (1994, page 691)

  f_sammler = rep(5,T-K)
  for(t in (K+1):T){
    f_sammler[t-K] = sum(zeta_t1_t[,t-1]*eta[,t]) #computes the value of the log-likelihood based on the formula in Hamilton (1994, page 692)
  }

  Log_Likelihood_value = sum(log(f_sammler + 0.000000001)) #add a small number to avoid zeros due to rounding or (early) extreme parameter values
  return(Log_Likelihood_value)
}
