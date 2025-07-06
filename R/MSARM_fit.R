
#' MSARM.fit
#' MSARM.fit allows its user to estimate Markov-Switching AR models, as they are described in Hamilton (1994, page 690-698). Thereby an arbitrary set of parameters
#' is allowed to switch. The estimation process strongly orients itself on the EM algorithm described in Hamilton (1990) and is described in more detail in Müller (2025, page 7-21).
#' @param Y_T Data (The full time series that is to be analyzed)
#' @param K Lag-order of the AR model
#' @param N Number of underlying Regimes
#' @param Switcher A vector of TRUE/FALSE indicating which parameters switch (intercept, coefficients, error term variance), the vector should have the length K+2
#' @param m Number of time periods the conditional likelihood function conditions on (generally set to 1)
#' @param threshold_value Threshold for assigning the estimated regimes (generally set to 0.5)
#' @param max_value Number of iterations for the EM algorithm (a higher iteration number increases the estimation quality)
#' @param R_value Number of random starting points for the EM algorithm (a higher number of random starting points increases estimation quality)
#' @param Crit_value Performance metric for choosing between the different optima. It is advised to run the function ones with "LV" and once with "RCM", if both estimation approaches lead to vastly different results a third run with Entropy can help to decide.
#' @param all.plot Whether plots from all optimisation attempts should be printed
#'
#' @return Returns a list object containing the estimation results regarding the Markov-Switching AR model. The output includes
#' the estimate of the P matrix, the estimate of the alpha matrix, the smoothed inference regarding the regines, the value of the log-likelihood function
#' the computed RSS value, the computed RCM value and the computed Entropy value. Furthermore the list contains the original time series, the in-sample fit and the residuals.
#'
#' @export
#'
#' @examples MSARM.fit(Y_T,K,N,Switcher) would be the standard implementation of MSARM.fit
MSARM.fit = function(Y_T,K,N,Switcher,m = 1,threshold_value = 0.5, max_value = 250, R_value = 5, Crit_value = "LV",all.plot = FALSE){


  #Based on the form of the switching vector it is decided whether Example 3, Example 4 or Example 5 is utilized as estimation framework
  if(sum(Switcher == c(rep(TRUE,K+1),FALSE)) == K+2){
    Example4_EM_Algorithm_function_lastbeta(Y_T = Y_T,
                                            K = K,
                                            N = N,
                                            m = m,
                                            threshold = threshold_value,
                                            max = max_value,
                                            R = R_value,
                                            Crit = Crit_value,
                                            all.plot = all.plot,
                                            Switcher = Switcher)
  }else{
    if(sum(Switcher == c(FALSE,rep(TRUE,K),FALSE))== K+2){
      Example4_EM_Algorithm_function_lastbeta(Y_T = Y_T,
                                              K = K,
                                              N = N,
                                              m = m,
                                              threshold = threshold_value,
                                              max = max_value,
                                              R = R_value,
                                              Crit = Crit_value,
                                              all.plot = all.plot,
                                              Switcher = Switcher)
    }else{
      if(sum(Switcher == c(TRUE,rep(FALSE,K),FALSE)) == K+2){
        Example4_EM_Algorithm_function_lastbeta(Y_T = Y_T,
                                                K = K,
                                                N = N,
                                                m = m,
                                                threshold = threshold_value,
                                                max = max_value,
                                                R = R_value,
                                                Crit = Crit_value,
                                                all.plot = all.plot,
                                                Switcher = Switcher)
      }else{
        if(sum(Switcher == c(rep(TRUE,K+2))) == K+2){
          Example3_EM_Algorithm_function_lastbeta(Y_T = Y_T,
                                                  K = K,
                                                  N = N,
                                                  m = m,
                                                  threshold = threshold_value,
                                                  max = max_value,
                                                  R = R_value,
                                                  Crit = Crit_value,
                                                  all.plot = all.plot,
                                                  Switcher = Switcher)
        }else{
          if(Switcher[K+2] == FALSE){
            Example4_EM_Algorithm_function_lastbeta(Y_T = Y_T,
                                                    K = K,
                                                    N = N,
                                                    m = m,
                                                    threshold = threshold_value,
                                                    max = max_value,
                                                    R = R_value,
                                                    Crit = Crit_value,
                                                    all.plot = all.plot,
                                                    Switcher = Switcher)
          }else{
            Example5_EM_Algorithm_function_lastbeta(Y_T = Y_T,
                                                    K = K,
                                                    N = N,
                                                    m = m,
                                                    threshold = threshold_value,
                                                    max = max_value,
                                                    R = R_value,
                                                    Crit = Crit_value,
                                                    all.plot = all.plot,
                                                    Switcher = Switcher)
          }
        }
      }
    }
  }
}
