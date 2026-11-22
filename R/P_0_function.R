#' @title Random Initial Transition Matrix
#' @param Y_T,K Retained for compatibility.
#' @param N Number of regimes.
#' @return A column-stochastic N by N transition matrix.
#' @export
P_0_function <- function(Y_T, N, K) {
  P <- matrix(0, N, N)
  for (i in seq_len(N)) {
    P[, i] <- stats::runif(N, 0.1, 0.25)
    P[i, i] <- stats::runif(1L, 0.75, 1)
  }
  sweep(P, 2L, colSums(P), "/")
}
