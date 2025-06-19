#' Title
#'
#' @param Y_T Data
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return A random starting value of zeta (regime probability) is generated
#' @export
#'
#' @examples zeta_1_0_function(Y_T,N,K)
zeta_1_0_function = function(Y_T,N,K){
  zeta = matrix(0,ncol = 1,nrow = N)

  p_col = runif(N,0,1)
  p_col_standardised = p_col/sum(p_col)
  zeta = p_col_standardised

  return(zeta)
}
