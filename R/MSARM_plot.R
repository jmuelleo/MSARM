#' MSARM.plot()
#'
#' @param erg_MSARM.predict output from MSARM.predict
#' @param conf if set to TRUE and MSARM.predict utilized boot = TRUE, then confidence intervals for the forecasts will be additionally plotted
#'
#' @return MSARM.plot allows its user to plot the forecasts from MSARM.predict to see how the time series is expected to behave in future. Furthermore
#' MSARM.plot allows its user if conf = TRUE and boot = TRUE (in MSARM.predict) to additionally plot the bootstrap confidence intervals for the forecasts
#' @export
#'
#' @examples MSARM.plot(erg_MSARM.predict,conf = TRUE)
MSARM.plot = function(erg_MSARM.predict,conf = TRUE){

  Y_T = erg_MSARM.predict$Y_T
  fcast_boot = erg_MSARM.predict$fcast_boot
  n.ahead = erg_MSARM.predict$n.ahead
  l = length(Y_T)
  total_series = ts(c(Y_T,fcast_boot[,1]))
  T = length(total_series)
  time_new = time(total_series)[(T-n.ahead):T]
  plot(ts(Y_T),lty = 1)
  time_new_poly = time(total_series)[(T-n.ahead+1):T]
  xx = c(time_new_poly,rev(time_new_poly))

  if(conf == TRUE){
    for(i in 1:((dim(fcast_boot)[2]-1)/2)){
      yy = c(fcast_boot[,i+1],rev(fcast_boot[,dim(fcast_boot)[2]-i+1]))
      polygon(xx,yy,col = gray(seq(0.8, 0.3, length.out = ((dim(fcast_boot)[2]-1)/2))[i], alpha = 0.8) )
    }
  }

  lines(time_new,c(Y_T[l],fcast_boot[,1]), col = "blue", lwd = 2)

  if(conf == TRUE){
    plot(ts(Y_T),lty = 1,xlim = c((T-n.ahead-50),T))
    time_new_poly = time(total_series)[(T-n.ahead+1):T]
    xx = c(time_new_poly,rev(time_new_poly))
    for(i in 1:((dim(fcast_boot)[2]-1)/2)){
      yy = c(fcast_boot[,i+1],rev(fcast_boot[,dim(fcast_boot)[2]-i+1]))
      polygon(xx,yy,col = gray(seq(0.8, 0.3, length.out = ((dim(fcast_boot)[2]-1)/2))[i], alpha = 0.8) )
    }
    lines(time_new,c(Y_T[l],fcast_boot[,1]), col = "blue", lwd = 2)
  }
}
