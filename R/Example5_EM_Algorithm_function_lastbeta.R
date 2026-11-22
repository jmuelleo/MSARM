

#' @title Example5_EM_Algorithm_function_lastbeta
#' @description Example5_EM_Algorithm_function_lastbeta allows its user to implement the EM algorithm for estimating the parameters in Example 5 setups, see Müller (2025, page 21-24)
#' @param Y_T Data
#' @param K Lag-order of the AR model
#' @param N Number of underlying regimes
#' @param m Transition counts start at observation K + m + 1. The default 1 includes every transition in the post-lag sample.
#' @param threshold Threshold for assigning the predicted regimes
#' @param max Number of iterations of the EM algorithm
#' @param R Number of random starting points for the EM algorithm
#' @param Crit Performance metric for selecting one of the optima
#' @param all.plot Whether plots from a optimisation attempts should be printed
#' @param Switcher Switching vector indicating for all parameters whether they switch or not
#'
#' @return Returns a list object containing the results from the EM algorithm estimation of the Markov-Switching model
#' @export
#'
Example5_EM_Algorithm_function_lastbeta <- function(
    Y_T, K, N, m, threshold = 0.5, max = 500, R = 15, Crit = "LV",
    all.plot = FALSE, Switcher) {
  msarm_fit(Y_T, K, N, Switcher, m, threshold, max, R, Crit, all.plot)
}
