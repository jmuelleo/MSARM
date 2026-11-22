
#' @title Log_Likelihood_function
#' @description Log_Likelihood_function allows its user to compute the value of the log-likelihood function based on Hamilton (1994, page 692)
#' @param alpha alpha vector from the previous iteration
#' @param P P matrix from the previous iteration
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return Returns the value of the log-likelihood function at a specific theta
#' @export
#'
Log_Likelihood_function <- function(alpha, P, Y_T, N, K) {
  msarm_filter(alpha, P, Y_T, N, K)$loglik
}
