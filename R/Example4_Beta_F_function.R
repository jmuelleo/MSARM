
#' @title Example4_Beta_F_function
#' @description Example4_Beta_F_function allows its user to estimate the fixed coefficients vector called beta_F, based on the results from the previous iteration and the formulas in Müller (2025, page 21)
#' @param Y_T Data
#' @param alpha alpha matrix from the previous iteration
#' @param P P matrix from the previous iteration
#' @param Switcher Switching vector indicating which parameters switch
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Retained for compatibility. This update uses the full post-lag sample.
#'
#' @return Returns a vector representing beta_F based on the results from the last iteration
#' @export
#'
Example4_Beta_F_function <- function(Y_T, alpha, P, Switcher, N, K, m) {
  fixed <- which(!Switcher[seq_len(K + 1L)])
  if (!length(fixed)) {
    return(0)
  }
  updated <- msarm_update_from_series(Y_T, alpha, P, Switcher, N, K)
  matrix(updated[1L, fixed], ncol = 1L)
}
