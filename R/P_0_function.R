#' Title
#'
#' @param Y_T Data
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return A random starting value of P is generated
#' @export
#'
#' @examples P_0_function(Y_T,N,K)
P_0_function = function(Y_T,N,K){
  T = length(Y_T)
  P = matrix(0,ncol = N,nrow =N)
  for(i in 1:N){
    p_col = runif(N,0.1,0.25)
    p_col[i] = runif(1,0.75,1)
    p_col_standardised = p_col/sum(p_col)
    P[,i] = p_col_standardised
  }
  return(P)
}
