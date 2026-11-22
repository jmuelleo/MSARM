
#' @title Example4_alpha_next_function_lastbeta
#' @description Example4_alpha_next_function_lastbeta allows its user to estimate the alpha matrix of the next iteration for Example 4 setups, based on Müller (2025, page 19-21)
#' @param Y_T Data
#' @param alpha alpha matrix from the last iteration
#' @param P P matrix from the last iteration
#' @param Switcher Switching vector which indicates which parameter switch
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Retained for compatibility. This update uses the full post-lag sample.
#'
#' @return Returns a matrix of the dimensions (N x K+2)
#' @export
#'
Example4_alpha_next_function_lastbeta <- function(
    Y_T, alpha, P, Switcher, N, K, m) {
  msarm_update_from_series(Y_T, alpha, P, Switcher, N, K)
}
