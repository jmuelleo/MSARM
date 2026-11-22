
#' @title alpha_0_function
#' @description A random starting value of alpha is generated so that the EM algorithm can start
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return Returns an N by (K + 2) matrix containing randomly generated values for (intercept, coefficients, error term variance)
#' @export
#'
#'
alpha_0_function <- function(Y_T, N, K) {
  spread <- stats::sd(Y_T)
  cbind(stats::runif(N, min(Y_T), max(Y_T)),
        matrix(stats::rnorm(N * K, sd = 0.25 / sqrt(K)), N, K),
        stats::runif(N, 0.25, 1) * spread^2)
}
