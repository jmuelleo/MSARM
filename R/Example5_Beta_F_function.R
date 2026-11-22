
#' @title Example5_Beta_F_function
#' @description Example5_Beta_F_function allows its user to compute beta_F as part of the EM algorithm for Example 5 setups, for details see Müller (2025, page 23-24)
#' @param Y_T Data
#' @param alpha alpha matrix from the previous iteration
#' @param P P matrix from the previous iteration
#' @param Switcher Switching vector indicating whether a parameter is fixed or switches
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Retained for compatibility. This update uses the full post-lag sample.
#'
#' @return Returns a vector representing the estimates for beta_F based on the results from the last iteration
#' @export
#'
Example5_Beta_F_function <- function(Y_T, alpha, P, Switcher, N, K, m) {
  fixed <- which(!Switcher[seq_len(K + 1L)])
  if (!length(fixed)) {
    return(0)
  }
  updated <- msarm_update_from_series(Y_T, alpha, P, Switcher, N, K)
  matrix(updated[1L, fixed], ncol = 1L)
}
