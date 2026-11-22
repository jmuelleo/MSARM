
#' @title P_ij_next_function
#' @description P_ij_next_function allows its user to estimate the P matrix of the next iteration as part of the EM algorithm, as described in Hamilton (1994, page 695)
#' @param alpha alpha matrix of the last iteration
#' @param P P matrix of the last iteration
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Transition counts start at observation K + m + 1. The default 1 includes every transition in the post-lag sample.
#'
#' @return Returns a matrix of the dimensions (N x N)
#' @export
#'
P_ij_next_function <- function(alpha, P, Y_T, N, K, m) {
  filter <- msarm_filter(alpha, P, Y_T, N, K)
  smoothed <- msarm_smooth(filter, P, K)
  msarm_transition(P, filter, smoothed, K, m)
}

msarm_transition <- function(P, filter, smoothed, K, m) {
  first <- K + m + 1L
  last <- ncol(smoothed)
  if (first > last) {
    return(P)
  }
  times <- seq.int(first, last)
  ratio <- msarm_probability_ratio(smoothed[, times, drop = FALSE],
                                   filter$t1t[, times - 1L, drop = FALSE])
  counts <- P * tcrossprod(ratio, filter$tt[, times - 1L, drop = FALSE])
  totals <- colSums(counts)
  occupied <- totals > 0
  P[, occupied] <- sweep(counts[, occupied, drop = FALSE], 2L,
                         totals[occupied], "/")
  P
}
