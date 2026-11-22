msarm_filter <- function(alpha, P, Y_T, N, K, design = msarm_design(Y_T, K)) {
  n <- length(Y_T)
  log_density <- msarm_log_density(alpha, Y_T, N, K, design)
  filtered <- predicted <- matrix(NA_real_, N, n)
  loglik <- numeric(n - K)
  predicted[, K] <- rep(1 / N, N)

  for (i in seq_along(loglik)) {
    time <- K + i
    log_weights <- log(predicted[, time - 1L]) + log_density[, i]
    offset <- max(log_weights)
    if (!is.finite(offset)) {
      stop("The conditional density is non-finite in every reachable regime.",
           call. = FALSE)
    }
    weights <- exp(log_weights - offset)
    total <- sum(weights)
    filtered[, time] <- weights / total
    predicted[, time] <- P %*% filtered[, time]
    loglik[i] <- offset + log(total)
  }
  list(tt = filtered, t1t = predicted, loglik = sum(loglik))
}

#' @title Filtered and Predicted Regime Probabilities
#' @rdname regime-filter
#' @description zeta_Yt_function allows its user to get optimal inference over the regimes as described in Hamilton (1994, page 692)
#' @param alpha alpha matrix from the last iteration
#' @param P P matrix from the last iteration
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#'
#' @return Matrices tt (filtered) and t1t (one-step predicted), with regimes in rows. Column t of t1t predicts time t + 1. Other presample entries are NA.
#' @export
#'
zeta_Yt_function <- function(alpha, P, Y_T, N, K) {
  result <- msarm_filter(alpha, P, Y_T, N, K)
  result[c("tt", "t1t")]
}
