
#' @title MSARM.fit
#' @description MSARM.fit allows its user to estimate Markov-Switching AR models, as they are described in Hamilton (1994, page 690-698). Thereby an arbitrary set of parameters
#' is allowed to switch. The estimation process strongly orients itself on the EM algorithm described in Hamilton (1990) and is described in more detail in Müller (2025, page 7-21).
#' @param Y_T Finite numeric vector or univariate time series.
#' @param K Lag-order of the AR model
#' @param N Number of underlying Regimes
#' @param Switcher A vector of TRUE/FALSE indicating which parameters switch (intercept, coefficients, error term variance), the vector should have the length K+2
#' @param m Transition counts start at observation K + m + 1. The default 1 includes every transition in the post-lag sample.
#' @param threshold_value Horizontal probability reference line when all.plot is TRUE.
#' @param max_value Number of iterations for the EM algorithm (a higher iteration number increases the estimation quality)
#' @param R_value Number of random starting points for the EM algorithm (a higher number of random starting points increases estimation quality)
#' @param Crit_value Criterion used to select a start. "LV" maximizes log-likelihood. "RSS", "RCM", and "Entropy" minimize the named statistic.
#' @param all.plot If TRUE, plot regime probabilities and fitted values for each start.
#'
#' @details Regime probabilities begin from a uniform distribution.
#'  Each start runs max_value iterations. Innovation variances are bounded
#'  below by 1e-10 times the sample variance.
#' @return A list with P (column-stochastic transitions from column to row),
#'  alpha (regime rows with intercept, AR coefficients, and variance),
#'  zeta_t_T (smoothed probabilities for the post-lag sample in columns),
#'  LV (conditional log-likelihood), RSS (residual sum of squares), RCM
#'  (regime classification measure on 0 to 100), Entropy (normalized Shannon
#'  entropy), Y_T (input), Yhat_insample (fitted values for the
#'  most probable smoothed regime), and residuals.
#' @references Baele, L. (2005). Volatility spillover effects in European
#'  equity markets. Journal of Financial and Quantitative Analysis, 40(2),
#'  373-401. doi:10.1017/S0022109000002350.
#' @examples
#' set.seed(12)
#' y <- arima.sim(list(ar = 0.5), n = 80)
#' fit <- MSARM.fit(y, 1, 2, c(TRUE, FALSE, TRUE),
#'                  max_value = 10, R_value = 1)
#' fit$P
#' @export
#'
MSARM.fit <- function(Y_T, K, N, Switcher, m = 1, threshold_value = 0.5,
                      max_value = 250, R_value = 5, Crit_value = "LV",
                      all.plot = FALSE) {
  msarm_fit(Y_T, K, N, Switcher, m, threshold_value, max_value, R_value,
            Crit_value, all.plot)
}
