#' MSARM.predict
#'
#' @param erg_MSARM.fit Result of an MSARM.fit estimation
#' @param n.ahead Number of periods to be forecasted ahead
#' @param boot Whether bootstrapping shall be utilized (allows for creating confidence intervals)
#' @param levels If bootstrapped, the levels for the confidence intervals
#' @param L If bootstrapped the number of bootstrap estimates
#'
#' @return MSARM.predict allows its user to predict an MSARM.fit output for n.ahead periods ahead. Furthermore confidence intervals can be created by setting
#' boot = TRUE, in that case L bootstrap estimates and intervals will be created.
#' @export
#'
#' @examples MSARM.predict(erg_MSARM.fit, n.ahead, boot = FALSE, levels = c(0.95,0.9,0.8,0.7,0.6), L = 10000)
MSARM.predict = function(erg_MSARM.fit, n.ahead, boot = FALSE, levels = c(0.95,0.9,0.8,0.7,0.6), L = 10000){

  zeta_t_T_1 = erg_MSARM.fit$zeta_t_T
  n = length(zeta_t_T_1[1,])
  zeta_T_T_1 = zeta_t_T_1[,n]

  P = erg_MSARM.fit$P


  Rhat = matrix(0,nrow = N, ncol = n.ahead)
  Rhat[,1] = P%*%zeta_T_T_1
  for(i in 2:n.ahead){
    Pm = P
    for (j in 2:i){
      Pm = Pm %*% P
    }
    Rhat[,i] = Pm %*%zeta_T_T_1
  }



  alpha = erg$alpha
  K = length(alpha[1,])-2
  coef_erg = alpha[,-(K+2)]
  Y_T = erg_MSARM.fit$Y_T

  fcast = rep(0,n.ahead)
  Y = c(Y_T)
  n = length(Y)
  h = rep(0,N)
  for(i in 1:N){
    h[i] = c(1,Y[n:(n-K+1)])%*%coef_erg[i,]
  }
  fcast[1] = t(h)%*%Rhat[,1]

  for(m in 2:n.ahead){
    Y = c(Y_T,fcast[1:(m-1)])
    n = length(Y)
    h = rep(0,N)
    for(i in 1:N){
      h[i] = c(1,Y[n:(n-K+1)])%*%coef_erg[i,]
    }
    fcast[m] = t(h)%*%Rhat[,m]
  }



  output = list(zeta_hat = Rhat,
                fcast = fcast,
                Y_T = Y_T,
                n.ahead = n.ahead)



  residuals = erg_MSARM.fit$residuals


  if(boot == TRUE){

    fcast_boot = matrix(0,nrow = n.ahead, ncol = L)
    for(l in 1:L){
      Y = c(Y_T)
      n = length(Y)
      h = rep(0,N)
      for(i in 1:N){
        h[i] = c(1,Y[n:(n-K+1)])%*%coef_erg[i,]
      }
      fcast_boot[1,l] = t(h)%*%Rhat[,1] + sample(residuals,1)

      for(m in 2:n.ahead){
        Y = c(Y_T,fcast[1:(m-1)])
        n = length(Y)
        h = rep(0,N)
        for(i in 1:N){
          h[i] = c(1,Y[n:(n-K+1)])%*%coef_erg[i,] + sample(residuals,1)
        }
        fcast_boot[m,l] = t(h)%*%Rhat[,m]
      }
    }

    quantile_vector = c((1-levels)/2,rev(levels+(1-levels)/2))
    erg_boot = matrix(0,nrow = n.ahead,ncol = (length(levels)*2)+1)
    colnames(erg_boot) = c("fcast_boot",quantile_vector)
    for(m in 1:n.ahead){

      erg_boot[m,1] = mean(fcast_boot[m,])

      for(i in 1:length(quantile_vector)){
        erg_boot[m,i+1] = quantile(fcast_boot[m,], probs = quantile_vector[i])
      }
    }

    output = list(zeta_hat = Rhat,
                  fcast = fcast,
                  fcast_boot = erg_boot,
                  Y_T = Y_T,
                  n.ahead = n.ahead)


  }

  return(output)

}
