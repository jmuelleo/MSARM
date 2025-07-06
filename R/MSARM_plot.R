
#' @title MSARM.plot
#' @description MSARM.plot creates forecast plots for MSARM.fit results
#' @param res_MSARM.predict An object saving the results from MSARM.predict
#' @param conf Set to TRUE if MSARM.predict was implemented with the option boot = TRUE, then bootstrap confidence intervals are plotted
#' @param start Start of the time series (format as in in the window() command), therefore (year,quarter) for quarterly data for example
#' @param freq Frequency of the data, i.e. 4 for quarterly, 12 for monthly data.
#'
#' @return Returns forecast plots
#' @export
#'
#' @examples MSARM.plot(res_MSARM.predict, start = c(1967,1), freq = 4) would give a standard forecast plot with the original data starting in 1967, first quarter
MSARM.plot = function(res_MSARM.predict,conf = FALSE, start = c(1,1), freq = 1){

  Y_T = res_MSARM.predict$Y_T #Original total time series
  fcast = res_MSARM.predict$fcast #forecast results from MSARM.predict
  n.ahead = res_MSARM.predict$n.ahead #number of forecasts
  l = length(Y_T)

  total_series = ts(c(Y_T,fcast), start = c(start), frequency = freq) #add forecasts to the original time series
  T = length(total_series)
  time_new = time(total_series)[(T-n.ahead):T]
  plot(ts(Y_T, start = start, freq = freq),lty = 1, xlim = c(range(time(total_series)))) #plot original time series
  time_new_poly = time(total_series)[(T-n.ahead+1):T]
  xx = c(time_new_poly,rev(time_new_poly))

  #add bootstrap forecast and bootstrap confidence intervals
  if(conf == TRUE){
    fcast_boot = res_MSARM.predict$fcast_boot
    for(i in 1:((dim(fcast_boot)[2]-1)/2)){
      yy = c(fcast_boot[,i+1],rev(fcast_boot[,dim(fcast_boot)[2]-i+1]))
      polygon(xx,yy,col = gray(seq(0.8, 0.3, length.out = ((dim(fcast_boot)[2]-1)/2))[i], alpha = 0.8) )
    }
  lines(time_new,c(Y_T[l],fcast_boot[,1]), col = "blue", lwd = 2)
  }else{
    lines(time_new,c(Y_T[l],fcast), col = "blue", lwd = 2) #if boot is set to FALSE plot the standard forecast
  }


 #plot only the last 30 observations and the bootstrap forecasts with confidence intervals
  if(conf == TRUE){
    fcast_boot = res_MSARM.predict$fcast_boot
    plot(ts(Y_T, start = start, freq = freq),lty = 1,xlim = c((time(total_series)[c((T-30-n.ahead),T)])))
    time_new_poly = time(total_series)[(T-n.ahead+1):T]
    xx = c(time_new_poly,rev(time_new_poly))
    for(i in 1:((dim(fcast_boot)[2]-1)/2)){
      yy = c(fcast_boot[,i+1],rev(fcast_boot[,dim(fcast_boot)[2]-i+1]))
      polygon(xx,yy,col = gray(seq(0.8, 0.3, length.out = ((dim(fcast_boot)[2]-1)/2))[i], alpha = 0.8) )
    }
    lines(time_new,c(Y_T[l],fcast_boot[,1]), col = "blue", lwd = 2)
  }
}
