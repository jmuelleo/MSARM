#' @title Draw Random Regime Weights
#' @description Generates normalized uniform draws. The fit and filtering
#'   routines use deterministic uniform initial probabilities.
#' @param Y_T,K Retained for compatibility.
#' @param N Number of regimes.
#' @return A probability vector of length N.
#' @export
zeta_1_0_function <- function(Y_T, N, K) {
  weights <- stats::runif(N)
  weights / sum(weights)
}
