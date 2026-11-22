
#' @title MSARM.predict
#' @description MSARM.predict allows its user to make forecasts utilizing an earlier estimated Markov-Switching AR model, the underlying theory is described in Müller (2025, page 24)
#' @param res_MSARM.fit An object containing the estimation results of MSARM.fit
#' @param n.ahead Number of time periods to forecast ahead
#' @param boot Whether bootstrap forecasts and bootstrap prediction intervals should be implemented instead of standard forecasts
#' @param levels The levels of the bootstrap prediction intervals, if boot is set to TRUE
#' @param L Number of bootstrap estimations computed
#'
#' @details Bootstrap paths simulate regime transitions and resample
#'  centered, standardized residuals. Forecast lags propagate along
#'  each path. Bounds are pointwise prediction intervals.
#' @return A list with zeta_hat, fcast, Y_T, and n.ahead. When boot
#'  is TRUE, fcast_boot contains the bootstrap mean and quantiles.
#' @export
#'
MSARM.predict <- function(res_MSARM.fit, n.ahead = 1, boot = FALSE,
                          levels = c(0.95, 0.9, 0.8, 0.7, 0.6), L = 10000) {
  n.ahead <- msarm_count(n.ahead, "n.ahead")
  if (!is.logical(boot) || length(boot) != 1L || is.na(boot)) {
    stop("boot must be TRUE or FALSE.", call. = FALSE)
  }
  P <- res_MSARM.fit$P
  alpha <- res_MSARM.fit$alpha
  Y_T <- res_MSARM.fit$Y_T
  probabilities <- res_MSARM.fit$zeta_t_T
  N <- nrow(P)
  K <- ncol(alpha) - 2L
  lags <- rev(utils::tail(as.numeric(Y_T), K))
  probability <- probabilities[, ncol(probabilities)]
  moments <- outer(probability, lags)
  coefficients <- alpha[, seq_len(K) + 1L, drop = FALSE]
  regime_forecast <- matrix(0, N, n.ahead)
  forecast <- numeric(n.ahead)

  # Each entry is E[Y at the given lag times I(current regime = j)].
  for (h in seq_len(n.ahead)) {
    probability <- drop(P %*% probability)
    predicted_moments <- P %*% moments
    current <- alpha[, 1L] * probability +
      rowSums(coefficients * predicted_moments)
    moments <- cbind(current, predicted_moments[, -K, drop = FALSE])
    regime_forecast[, h] <- probability
    forecast[h] <- sum(current)
  }
  output <- list(zeta_hat = regime_forecast, fcast = forecast, Y_T = Y_T,
                 n.ahead = n.ahead)
  if (!boot) {
    return(output)
  }

  L <- msarm_count(L, "L")
  if (!is.numeric(levels) || !length(levels) || any(!is.finite(levels)) ||
      any(levels <= 0 | levels >= 1)) {
    stop("levels must contain probabilities strictly between zero and one.",
         call. = FALSE)
  }
  levels <- sort(levels, decreasing = TRUE)
  bounds <- c((1 - levels) / 2, rev((1 + levels) / 2))
  sd <- sqrt(alpha[, K + 2L])
  regimes <- max.col(t(probabilities), ties.method = "first")
  errors <- as.numeric(res_MSARM.fit$residuals) / sd[regimes]
  errors <- errors - mean(errors)
  states <- sample.int(N, L, replace = TRUE,
                       prob = probabilities[, ncol(probabilities)])
  history <- matrix(rep(lags, each = L), L, K)
  intervals <- matrix(0, n.ahead, 1L + length(bounds),
                      dimnames = list(NULL, c("fcast_boot", bounds)))
  for (h in seq_len(n.ahead)) {
    for (j in seq_len(N)) {
      paths <- which(states == j)
      states[paths] <- -sample.int(N, length(paths), replace = TRUE,
                                   prob = P[, j])
    }
    states <- -states
    noise <- errors[sample.int(length(errors), L, replace = TRUE)]
    values <- alpha[states, 1L] +
      rowSums(coefficients[states, , drop = FALSE] * history) +
      sd[states] * noise
    history <- cbind(values, history[, -K, drop = FALSE])
    intervals[h, ] <- c(mean(values), stats::quantile(values, bounds))
  }
  output$fcast_boot <- intervals
  output[c("zeta_hat", "fcast", "fcast_boot", "Y_T", "n.ahead")]
}
