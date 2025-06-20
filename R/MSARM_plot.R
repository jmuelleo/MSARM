#' MSARM.plot()
#'
#' @param res_MSARM.predict output from MSARM.predict
#' @param conf if set to TRUE and MSARM.predict utilized boot = TRUE, then confidence intervals for the forecasts will be additionally plotted
#'
#' @return MSARM.plot allows its user to plot the forecasts from MSARM.predict to see how the time series is expected to behave in future. Furthermore
#' MSARM.plot allows its user if conf = TRUE and boot = TRUE (in MSARM.predict) to additionally plot the bootstrap confidence intervals for the forecasts
#' @export
#'
#' @examples MSARM.plot(res_MSARM.predict,conf = TRUE)
MSARM.plot = function(res_MSARM.predict,conf = FALSE, start = c(1,1), freq = 1){

  Y_T = res_MSARM.predict$Y_T
  fcast = res_MSARM.predict$fcast
  n.ahead = res_MSARM.predict$n.ahead
  l = length(Y_T)

  total_series = ts(c(Y_T,fcast), start = c(start), frequency = freq)
  T = length(total_series)
  time_new = time(total_series)[(T-n.ahead):T]
  plot(ts(Y_T, start = start, freq = freq),lty = 1, xlim = c(range(time(total_series))))
  time_new_poly = time(total_series)[(T-n.ahead+1):T]
  xx = c(time_new_poly,rev(time_new_poly))

  if(conf == TRUE){
    fcast_boot = res_MSARM.predict$fcast_boot
    for(i in 1:((dim(fcast_boot)[2]-1)/2)){
      yy = c(fcast_boot[,i+1],rev(fcast_boot[,dim(fcast_boot)[2]-i+1]))
      polygon(xx,yy,col = gray(seq(0.8, 0.3, length.out = ((dim(fcast_boot)[2]-1)/2))[i], alpha = 0.8) )
    }
  lines(time_new,c(Y_T[l],fcast_boot[,1]), col = "blue", lwd = 2)
  }else{
    lines(time_new,c(Y_T[l],fcast), col = "blue", lwd = 2)
  }



  if(conf == TRUE){
    fcast_boot = res_MSARM.predict$fcast_boot
    plot(ts(Y_T, start = start, freq = freq),lty = 1,xlim = c((time(total_series)[(T-30-n.ahead),T])))
    time_new_poly = time(total_series)[(T-n.ahead+1):T]
    xx = c(time_new_poly,rev(time_new_poly))
    for(i in 1:((dim(fcast_boot)[2]-1)/2)){
      yy = c(fcast_boot[,i+1],rev(fcast_boot[,dim(fcast_boot)[2]-i+1]))
      polygon(xx,yy,col = gray(seq(0.8, 0.3, length.out = ((dim(fcast_boot)[2]-1)/2))[i], alpha = 0.8) )
    }
    lines(time_new,c(Y_T[l],fcast_boot[,1]), col = "blue", lwd = 2)
  }
}
