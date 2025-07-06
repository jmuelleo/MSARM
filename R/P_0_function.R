
#' P_0_function
#' P_0_function creates P_0 as part of a random starting point for the EM algorithm
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return
#' @export
#'
#' @examples P_0_function(Y_T,N,K)
P_0_function = function(Y_T,N,K){
  T = length(Y_T)
  P = matrix(0,ncol = N,nrow =N)
  #Random P_0 generation
  for(i in 1:N){
    p_col = runif(N,0.1,0.25)
    p_col[i] = runif(1,0.75,1)
    p_col_standardised = p_col/sum(p_col)
    P[,i] = p_col_standardised
  }
  return(P)
}
