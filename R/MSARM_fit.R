#' MSARM.fit
#'
#' @param Y_T Data
#' @param K Lag order of the AR-Model
#' @param N Number of underlying Regimes
#' @param Switcher Switching-vector, a vector of (TRUE/FALSE) of the length K+2, which indicates which parameters are supposed to switch
#' @param m Number of observations to condition on (generally set to 1)
#' @param threshold_value Threshold for assigning the underlying regime based on the estimated probability
#' @param max_value Number of iterations to go through
#' @param R_value Number of random starting points for the optimisation
#' @param Crit_value Performance Metric to use for choosing the optimisation result ("LV","RSS","RCM","Entropy")
#'
#' @return MSARM.fit allows its user to estimate Markov-Switching AR-Models.
#' The potential performance Metrics that could be used are:
#' "LV" -> Utilizes the estimated conditonal log-likelhood function for optima selection
#' "RSS" -> Utilizes the quality of the in-sample fit for optima selection
#' "RCM" -> Utilizes an Gini-Coefficient approach for optima selection
#' "Entropy" -> Utilizes an Entropy approach for optima selection
#' MSARM utilizes the EM-Algorithm outlined in Hamilton(1994) and Hamilton(1990).
#' Furthermore more explanation regarding the theory implemented and how it builds beyond Hamiltons work can be found in From Theory to Application: Developing MSARM, an R Package for Markov-Switching Autoregressive Models
#'
#' @export
#'
#' @examples MSARM.fit(Y_T,K,N,Switcher,m = 1,threshold_value = 0.5, max_value = 250, R_value = 5, Crit_value = "LV")
MSARM.fit = function(Y_T,K,N,Switcher,m = 1,threshold_value = 0.5, max_value = 250, R_value = 5, Crit_value = "LV",all.plot = FALSE){

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
