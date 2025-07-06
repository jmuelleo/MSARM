
#' MSARM.predict
#' MSARM.predict allows its user to make forecasts utilizing an earlier estimated Markov-Switching AR model, the underlying theory is described in Müller (2025, page 24)
#' @param res_MSARM.fit An object containing the estimation results of MSARM.fit
#' @param n.ahead Number of time periods to forecast ahead
#' @param boot Whether bootstrap forecasts and bootstrap confidence intervals should be implemented instead of standard forecasts
#' @param levels The levels of the boostrap confidence intervals, if boot is set to TRUE
#' @param L Number of boostrap estimations computed
#'
#' @return Returns a list object containg the results of the forecasts based on the estimated Markov-Switching model. The list object contains
#' the predicted future regime chance, the forecasts, the original time series and the number of time periods forecasted ahead. If boot is set to true, then also bootstrap confidence intervals are given.
#' @export
#'
#' @examples MSARM.predict(res_MSARM.fit, n.ahead = 4) would give for a quarterly time series the forecast for the next year.
MSARM.predict = function(res_MSARM.fit, n.ahead = 1, boot = FALSE, levels = c(0.95,0.9,0.8,0.7,0.6), L = 10000){

  N = dim(res_MSARM.fit$P)[1]
  zeta_t_T_1 = res_MSARM.fit$zeta_t_T #smoothed inference over the regimes
  n = length(zeta_t_T_1[1,])
  zeta_T_T_1 = zeta_t_T_1[,n]

  P = res_MSARM.fit$P #P matrix


  #create the m-step ahead predicted regime chance
  Rhat = matrix(0,nrow = N, ncol = n.ahead)
  Rhat[,1] = P%*%zeta_T_T_1
  if(n.ahead > 1){
  for(i in 2:n.ahead){
    Pm = P
    for (j in 2:i){
      Pm = Pm %*% P
    }
    Rhat[,i] = Pm %*%zeta_T_T_1
  }
}


  alpha = res_MSARM.fit$alpha #alpha matrix
  K = length(alpha[1,])-2
  coef_res = alpha[,-(K+2)] #estimated coefficients
  Y_T = res_MSARM.fit$Y_T #original time series

  #Create one-step ahead forecast
  fcast = rep(0,n.ahead)
  Y = c(Y_T)
  n = length(Y)
  h = rep(0,N)
  for(i in 1:N){
    h[i] = c(1,Y[n:(n-K+1)])%*%coef_res[i,]
  }
  fcast[1] = t(h)%*%Rhat[,1]

  #If more than one-step ahead forecasts are implemented:
  if(n.ahead > 1){
  for(m in 2:n.ahead){
    Y = c(Y_T,fcast[1:(m-1)])
    n = length(Y)
    h = rep(0,N)
    for(i in 1:N){
      h[i] = c(1,Y[n:(n-K+1)])%*%coef_res[i,] #m-step ahead conditional forecast
    }
    fcast[m] = t(h)%*%Rhat[,m] #create m-step ahead forecast
  }
}


  #save results
  output = list(zeta_hat = Rhat,
                fcast = fcast,
                Y_T = Y_T,
                n.ahead = n.ahead)



  residuals = res_MSARM.fit$residuals #in-sample residuals


  #run bootstrap algorithm as described in Müller (2025, page 28) if boot is set to TRUE
  if(boot == TRUE){

    fcast_boot = matrix(0,nrow = n.ahead, ncol = L)
    for(l in 1:L){
      Y = c(Y_T)
      n = length(Y)
      h = rep(0,N)
      for(i in 1:N){
        h[i] = c(1,Y[n:(n-K+1)])%*%coef_res[i,]
      }
      fcast_boot[1,l] = t(h)%*%Rhat[,1] + sample(residuals,1) #one-step ahead forecast

      if(n.ahead > 1){
      for(m in 2:n.ahead){
        Y = c(Y_T,fcast[1:(m-1)])
        n = length(Y)
        h = rep(0,N)
        for(i in 1:N){
          h[i] = c(1,Y[n:(n-K+1)])%*%coef_res[i,] + sample(residuals,1)
        }
        fcast_boot[m,l] = t(h)%*%Rhat[,m] #m-step ahead forecasts
      }}
    }

    #create a suitable results matrix
    quantile_vector = c((1-levels)/2,rev(levels+(1-levels)/2))
    res_boot = matrix(0,nrow = n.ahead,ncol = (length(levels)*2)+1)
    colnames(res_boot) = c("fcast_boot",quantile_vector)
    for(m in 1:n.ahead){

      res_boot[m,1] = mean(fcast_boot[m,])

      for(i in 1:length(quantile_vector)){
        res_boot[m,i+1] = quantile(fcast_boot[m,], probs = quantile_vector[i])
      }
    }

    #save results
    output = list(zeta_hat = Rhat,
                  fcast = fcast,
                  fcast_boot = res_boot,
                  Y_T = Y_T,
                  n.ahead = n.ahead)


  }

  return(output)

}
