msarm_log_density <- function(alpha, Y_T, N, K, design = msarm_design(Y_T, K)) {
  means <- alpha[, seq_len(K + 1L), drop = FALSE] %*% t(design$X)
  matrix(stats::dnorm(rep(design$y, each = N), means,
                      sqrt(alpha[, K + 2L]), log = TRUE), nrow = N)
}

#' @title eta_function
#' @description The eta-vector, defined as in Hamilton (1994, page 691) is computed
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return An N by length(Y_T) matrix of conditional densities. Presample columns are NA.
#' @export
#'
eta_function <- function(alpha, Y_T, N, K) {
  cbind(matrix(NA_real_, N, K), exp(msarm_log_density(alpha, Y_T, N, K)))
}
