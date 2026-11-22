
#' @title Example4_EM_Algorithm_function_lastbeta
#' @description Example4_EM_Algorithm_function_lastbeta allows its user to run the EM algorithm for estimating the parameters of a Markov-Switching AR model of the form of Example 4, see Müller (2025, page 19-21)
#' @param Y_T Data
#' @param K Lag-order of the AR model
#' @param N Number of underlying regimes
#' @param m Transition counts start at observation K + m + 1. The default 1 includes every transition in the post-lag sample.
#' @param threshold Threshold for assigning the regimes
#' @param max Number of iterations of the EM algorithm
#' @param R Number of random starting points of the EM algorithm
#' @param Crit Performance metric used for choosing between the optima
#' @param all.plot Whether the plots from all optimisation attempts should be printed
#' @param Switcher Switching vector indicating for all parameters whether they switch or not
#'
#' @return Returns a list object containing the results from the EM algorithm estimation of the Markov-Switching AR model
#' @export
#'
Example4_EM_Algorithm_function_lastbeta <- function(
    Y_T, K, N, m, threshold = 0.5, max = 500, R = 15, Crit = "LV",
    all.plot = FALSE, Switcher) {
  msarm_fit(Y_T, K, N, Switcher, m, threshold, max, R, Crit, all.plot)
}
