msarm_probability_ratio <- function(smoothed, predicted) {
  ratio <- smoothed
  positive <- predicted > 0
  ratio[positive] <- smoothed[positive] / predicted[positive]
  ratio[!positive] <- 0
  ratio
}

msarm_smooth <- function(filter, P, K) {
  smoothed <- filter$tt
  n <- ncol(smoothed)
  if (n > K + 1L) {
    for (time in seq.int(n - 1L, K + 1L)) {
      ratio <- msarm_probability_ratio(smoothed[, time + 1L],
                                       filter$t1t[, time])
      weights <- filter$tt[, time] * drop(crossprod(P, ratio))
      smoothed[, time] <- weights / sum(weights)
    }
  }
  smoothed
}

#' @title Smoothed Regime Probabilities
#' @description zeta_YT_function allows its user to get smoothed inference over the regimes as described in Hamilton (1994, page 694)
#' @param alpha alpha matrix from the last iteration
#' @param P P matrix from the last iteration
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return An N by length(Y_T) matrix of smoothed probabilities. Presample columns are NA.
#' @export
#'
zeta_YT_function <- function(alpha, P, Y_T, N, K) {
  msarm_smooth(msarm_filter(alpha, P, Y_T, N, K), P, K)
}
