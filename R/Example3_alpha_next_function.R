
#' @title Example3_alpha_next_function
#' @description The next alpha for Example 3 setups (all parameters switch) is computed as part of the EM algorithm
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param P transpose of the transition matrix (of the Markov-Chain) from the last iteration
#' @param Y_T Data
#' @param N Number of regimes
#' @param K Lag-order of the AR model
#' @param m Retained for compatibility. This update uses the full post-lag sample.
#'
#' @return returns a matrix of the dimensions (N x K+2)
#' @export
#'
Example3_alpha_next_function <- function(alpha, P, Y_T, N, K, m) {
  msarm_update_from_series(Y_T, alpha, P, rep(TRUE, K + 2L), N, K)
}
