
#' @title Example5_alpha_next_function_lastbeta
#' @description Example5_alpha_next_function_lastbeta allows its user to estimate the alpha matrix of the next iteration as part of the EM algorithm for setups of the form of Example 5. A more detailed explanation can be found in Müller (2025, page 21-24)
#' @param Y_T Data
#' @param alpha alpha matrix from the previous iteration
#' @param P P matrix from the previous iteration
#' @param Switcher Switching vector indicating which parameters switch
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Retained for compatibility. This update uses the full post-lag sample.
#'
#' @return Returns a matrix of the dimensions (N x K+2)
#' @export
#'
Example5_alpha_next_function_lastbeta <- function(
    Y_T, alpha, P, Switcher, N, K, m) {
  msarm_update_from_series(Y_T, alpha, P, Switcher, N, K)
}
