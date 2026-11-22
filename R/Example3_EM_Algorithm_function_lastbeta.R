
#' @title Example3_EM_Algorithm_function_lastbeta
#' @description Example3_EM_Algorithm_function_lastbeta allows its user to apply the EM algorithm for estimating Markov-Switching models of the form of Example 3 (All Switching), as described in Müller (2025, page 18-19)
#' @param Y_T Data
#' @param K Lag-order of the AR model
#' @param N Number of underlying regimes
#' @param m Transition counts start at observation K + m + 1. The default 1 includes every transition in the post-lag sample.
#' @param threshold #threshold for assigning the predicted regimes
#' @param max #number of iterations for the EM algorithm
#' @param R #number of random starting points for the EM algorithm
#' @param Crit #performance metric for choosing one of the R optima
#' @param all.plot #whether plots from all R optimisation attempts should be printed
#' @param Switcher #Switching vector indicating which parameters switch
#'
#' @return returns a list object containg the results from the EM algorithm estimation of the Markov-Switching AR model
#' @export
#'
Example3_EM_Algorithm_function_lastbeta <- function(
    Y_T, K, N, m, threshold = 0.5, max = 500, R = 15, Crit = "LV",
    all.plot = FALSE, Switcher) {
  msarm_fit(Y_T, K, N, Switcher, m, threshold, max, R, Crit, all.plot)
}
