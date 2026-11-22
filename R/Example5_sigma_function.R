
#' @title Example5_sigma_function
#' @description Example5_sigma_function allows its user to esimate the error term standard deviation for Example 5 setups, based on the results from previous iterations and the formula presented in Müller (2025, page 24)
#' @param Y_T Data
#' @param alpha alpha matrix from the previous iteration
#' @param P P matrix from the previous iteration
#' @param Switcher Switching vector indicating whether a parameter switches or not
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Retained for compatibility. This update uses the full post-lag sample.
#'
#' @return Returns a vector representing the estimates for the standard deviation of the error term, based on the results from the last iteration
#' @export
#'
Example5_sigma_function <- function(Y_T, alpha, P, Switcher, N, K, m) {
  design <- msarm_design(Y_T, K)
  probabilities <- zeta_YT_function(alpha, P, Y_T, N, K)[, -seq_len(K),
                                                         drop = FALSE]
  residuals <- sweep(design$X %*% t(alpha[, seq_len(K + 1L), drop = FALSE]),
                     1L, design$y, "-")
  mass <- rowSums(probabilities)
  sigma2 <- alpha[, K + 2L]
  active <- mass > 0
  sigma2[active] <- rowSums(probabilities * t(residuals^2))[active] /
    mass[active]
  sqrt(pmax(sigma2, stats::var(Y_T) * 1e-10))
}
