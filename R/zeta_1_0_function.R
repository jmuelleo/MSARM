
#' zeta_1_0_function
#' zeta_1_0_function creates zeta_1_0 as part of the random starting point for the EM algorithm
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return
#' @export
#'
#' @examples zeta_1_0_function(Y_T,N,K)
zeta_1_0_function = function(Y_T,N,K){
  T = length(Y_T)
  zeta = matrix(0,ncol = 1,nrow = N)

  p_col = runif(N,0,1)
  p_col_standardised = p_col/sum(p_col) #random zeta_1_0 to start the optimal inference over the regimes, which is essential for the EM algorithm
  zeta = p_col_standardised

  return(zeta)
}
