#' Title
#'
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param Y_T Data
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return The eta-vector, defined as in Hamilton(1994) p.691 is computed
#' @export
#'
#' @examples eta_function(alpha,Y_T,N,K)
eta_function = function(alpha,Y_T,N,K){
  T = length(Y_T)
  eta = matrix(5,ncol = T,nrow = N)
  for(t in (K+1):T){ #Vorstichprobe, daher + K
    for(h in 1:N){
      eta[h,t] = dnorm(Y_T[t],
                       mean = c(1,Y_T[(t-1):(t-K)])%*%alpha[h,-(K+2)],
                       sd = sqrt(alpha[h,(K+2)]))

    }
  }
  return(eta + 0.000000001) #Die ersten Iterationen können sehr weit weg sein, da ein zufälliges alpha_0 gewählt wird, daher eine kleine Konstante, damit durch Rundung da keine echte 0 steht!
}
