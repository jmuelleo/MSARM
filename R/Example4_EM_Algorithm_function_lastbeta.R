#' Title
#'
#' @param Y_T Data
#' @param K Lag order of the AR-Model
#' @param N Number of underlying Regimes
#' @param m Number of observations to condition on
#' @param threshold Threshold for assigning the underlying regime based on the estimated probability
#' @param max Number of iterations to go through
#' @param R Number of random starting points
#' @param Crit Performance Metric to use for choosing the optimisation result
#'
#' @return This function brings, for Example 4 setups (Subset of the coefficient switches, error term variance is fixed), everything together and runs the EM-Algorithm
#' @export
#'
#' @examples Example4_EM_Algorithm_function_lastbeta(Y_T,K,N,m,threshold = 0.5,max = 500,R = 15, Crit = "LV")
Example4_EM_Algorithm_function_lastbeta = function(Y_T,K,N,m,threshold = 0.5,max = 500,R = 15, Crit = "LV", all.plot = FALSE){

  Output_list = as.list(1:R)

  LV_list = as.list(1:R)
  RSS_list = as.list(1:R)
  RCM_list = as.list(1:R)
  Entropy_list = as.list(1:R)
  for(r in 1:R){


    zetahat_1_0 = zeta_1_0_function(Y_T,N,K)
    alpha_0 = alpha_0_function(Y_T,N,K)
    P_0 = P_0_function(Y_T,N,K)
    T = length(Y_T)
    counter = 1
    P_Sammler = as.list(1:(max+1))
    alpha_Sammler = as.list(1:(max+1))
    P_Sammler[[1]] = P_0
    alpha_Sammler[[1]] = alpha_0
    while(counter <= max){

      P_Sammler[[counter+1]] = P_ij_next_function(alpha = alpha_Sammler[[counter]],P = P_Sammler[[counter]],Y_T = Y_T,N = N,K = K,m)
      alpha_Sammler[[counter+1]] = Example4_alpha_next_function_lastbeta(Y_T,alpha_Sammler[[counter]],P_Sammler[[counter]],Switcher,N,K,m)
      Log_value_l = Log_Likelihood_function(alpha_Sammler[[ifelse(counter-1>=1,counter-1,1)]],P_Sammler[[ifelse(counter-1>=1,counter-1,  1)]],Y_T,N,K)
      Log_value_n = Log_Likelihood_function(alpha_Sammler[[counter]],P_Sammler[[counter]],Y_T,N,K)

      Delta = abs(Log_value_n - Log_value_l)

      # if(Delta <= 0.00001){
      #   max = counter - 1
      #   break
      #}

      counter = counter + 1
    }
    Log_value = Log_Likelihood_function(alpha_Sammler[[counter]],P_Sammler[[counter]],Y_T,N,K)




    zetaout = zeta_YT_function(alpha_Sammler[[max+1]],P_Sammler[[max+1]],Y_T,N,K)[,-(1:(K))]
    zetaout_ts = ts(zetaout[1,]) #1 oder 2 wählen je nachdem, wie das Modell die Regime zuordnet
if(all.plot == TRUE){
    time_vals = time(zetaout_ts)
    plot(zetaout_ts, type = "l", col = "blue", lwd = 2,
         main = "Probability of Regime 2",
         xlab = "Year", ylab = "Probability")
    above_threshold = zetaout_ts > threshold
    for (i in which(above_threshold)) {
      rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
           time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
           col = rgb(1, 0, 0, alpha = 0.2), border = NA)
    }
    lines(zetaout_ts, col = "blue", lwd = 2)
    abline(h = threshold, col = "black", lty = 2)
}
    years = floor(time(zetaout_ts))

if(all.plot == TRUE){
    zetaout_ts_M = ts(t(zetaout))
    plot(zetaout_ts_M,plot.type = "single", col = c(1:N), lwd = 2)
    plot(zetaout_ts_M, col = c(1:N), lwd = 2)
}
    Predicted_Regime = apply(zetaout_ts_M,1,function(row){which(row == max(row))})

    Yhat_insample = rep(0,T-K)
    for(t in 1:(T-K)){
      for(j in 1:N){

        if(j == Predicted_Regime[t]){

          Yhat_insample[t] = c(c(1,Y_T[(t+K-1):(t)])%*%alpha_Sammler[[max+1]][j,-(K+2)])

        }

      }
    }
if(all.plot == TRUE){
    xplot = time(Y_T[-(1:K)])
    xxplot = c(xplot,rev(xplot))
    yyplot = c(Y_T[-(1:K)],rev(ts(Yhat_insample)))
    plot(cbind(Y_T[-(1:K)],ts(Yhat_insample)),plot.type = "single", col = c("black","blue"),lwd = c(1,2), lty = 1,main = "In Sample Fit", ylab = "Values")
    polygon(xxplot,yyplot,col = "lightblue",border = FALSE)
    for (i in which(above_threshold)) {
      rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
           time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
           col = rgb(1, 0, 0, alpha = 0.2), border = NA)
    }
    lines(ts(Yhat_insample),lwd = 2, col = "blue")
    lines(Y_T[-(1:K)],lwd = 1, col = "black")
    legend("bottomleft",fill = c("black","blue","red"),legend = c("Actual Time Series","In Sample Fit","Predicted Regime 2"))

}
    RSS = sum((Y_T[-(1:K)] - ts(Yhat_insample))^2)
    residuals = Y_T[-(1:K)] - Yhat_insample

    RCM = 100*(N^2)*mean(apply(zetaout,2,prod))


    Persistency = prod(diag(P_Sammler[[max+1]]))

    Entropy = 100*(N^2)*(-1)*mean(apply(zetaout,2,function(col){
      phat = col[1]
      return(phat*log(phat))
    }))

    Entropy = Entropy

    output = list(P = P_Sammler[[max+1]],
                  alpha = alpha_Sammler[[max+1]],
                  zeta_t_T = zetaout,
                  LV  = Log_value,
                  RSS = RSS,
                  RCM = RCM,
                  Entropy = Entropy,
                  Y_T = Y_T,
                  Yhat_insample = Yhat_insample,
                  residuals = residuals)


    print(alpha_Sammler[[max+1]])
    print(P_Sammler[[max+1]])

    Output_list[[r]] = output
    LV_list[[r]] = output$LV
    RSS_list[[r]] = output$RSS
    RCM_list[[r]] = output$RCM
    Entropy_list[[r]] = output$Entropy
  }
  if(Crit == "LV"){
    output = Output_list[[which.max(LV_list)]]
    P = output$P
    alpha = output$alpha
    zeta_out = output$zeta_t_T
    Y_T = output$Y_T
    Yhat_insample = output$Yhat_insample


    zetaout_ts = ts(zetaout[1,])
    time_vals = time(zetaout_ts)
    plot(zetaout_ts, type = "l", col = "blue", lwd = 2,
         main = "Probability of Regime 2",
         xlab = "Year", ylab = "Probability")
    above_threshold = zetaout_ts > threshold
    for (i in which(above_threshold)) {
      rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
           time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
           col = rgb(1, 0, 0, alpha = 0.2), border = NA)
    }
    lines(zetaout_ts, col = "blue", lwd = 2)
    abline(h = threshold, col = "black", lty = 2)


    zetaout_ts_M = ts(t(zetaout))
    plot(zetaout_ts_M,plot.type = "single", col = c(1:N), lwd = 2)
    plot(zetaout_ts_M, col = c(1:N), lwd = 2)





    xplot = time(Y_T[-(1:K)])
    xxplot = c(xplot,rev(xplot))
    yyplot = c(Y_T[-(1:K)],rev(ts(Yhat_insample)))
    plot(cbind(Y_T[-(1:K)],ts(Yhat_insample)),plot.type = "single", col = c("black","blue"),lwd = c(1,2), lty = 1,main = "In Sample Fit", ylab = "Values")
    polygon(xxplot,yyplot,col = "lightblue",border = FALSE)
    for (i in which(above_threshold)) {
      rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
           time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
           col = rgb(1, 0, 0, alpha = 0.2), border = NA)
    }
    lines(ts(Yhat_insample),lwd = 2, col = "blue")
    lines(Y_T[-(1:K)],lwd = 1, col = "black")
    legend("bottomleft",fill = c("black","blue","red"),legend = c("Actual Time Series","In Sample Fit","Predicted Regime 2"))
  }else{
    if(Crit == "RSS"){
      output = Output_list[[which.min(RSS_list)]]
      P = output$P
      alpha = output$alpha
      zeta_out = output$zeta_t_T
      Y_T = output$Y_T
      Yhat_insample = output$Yhat_insample


      zetaout_ts = ts(zetaout[1,])
      time_vals = time(zetaout_ts)
      plot(zetaout_ts, type = "l", col = "blue", lwd = 2,
           main = "Probability of Regime 2",
           xlab = "Year", ylab = "Probability")
      above_threshold = zetaout_ts > threshold
      for (i in which(above_threshold)) {
        rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
             time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
             col = rgb(1, 0, 0, alpha = 0.2), border = NA)
      }
      lines(zetaout_ts, col = "blue", lwd = 2)
      abline(h = threshold, col = "black", lty = 2)


      zetaout_ts_M = ts(t(zetaout))
      plot(zetaout_ts_M,plot.type = "single", col = c(1:N), lwd = 2)
      plot(zetaout_ts_M, col = c(1:N), lwd = 2)





      xplot = time(Y_T[-(1:K)])
      xxplot = c(xplot,rev(xplot))
      yyplot = c(Y_T[-(1:K)],rev(ts(Yhat_insample)))
      plot(cbind(Y_T[-(1:K)],ts(Yhat_insample)),plot.type = "single", col = c("black","blue"),lwd = c(1,2), lty = 1,main = "In Sample Fit", ylab = "Values")
      polygon(xxplot,yyplot,col = "lightblue",border = FALSE)
      for (i in which(above_threshold)) {
        rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
             time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
             col = rgb(1, 0, 0, alpha = 0.2), border = NA)
      }
      lines(ts(Yhat_insample),lwd = 2, col = "blue")
      lines(Y_T[-(1:K)],lwd = 1, col = "black")
      legend("bottomleft",fill = c("black","blue","red"),legend = c("Actual Time Series","In Sample Fit","Predicted Regime 2"))
    }else{
      if(Crit == "RCM"){
        output = Output_list[[which.min(RCM_list)]]
        P = output$P
        alpha = output$alpha
        zeta_out = output$zeta_t_T
        Y_T = output$Y_T
        Yhat_insample = output$Yhat_insample


        zetaout_ts = ts(zetaout[1,])
        time_vals = time(zetaout_ts)
        plot(zetaout_ts, type = "l", col = "blue", lwd = 2,
             main = "Probability of Regime 2",
             xlab = "Year", ylab = "Probability")
        above_threshold = zetaout_ts > threshold
        for (i in which(above_threshold)) {
          rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
               time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
               col = rgb(1, 0, 0, alpha = 0.2), border = NA)
        }
        lines(zetaout_ts, col = "blue", lwd = 2)
        abline(h = threshold, col = "black", lty = 2)


        zetaout_ts_M = ts(t(zetaout))
        plot(zetaout_ts_M,plot.type = "single", col = c(1:N), lwd = 2)
        plot(zetaout_ts_M, col = c(1:N), lwd = 2)





        xplot = time(Y_T[-(1:K)])
        xxplot = c(xplot,rev(xplot))
        yyplot = c(Y_T[-(1:K)],rev(ts(Yhat_insample)))
        plot(cbind(Y_T[-(1:K)],ts(Yhat_insample)),plot.type = "single", col = c("black","blue"),lwd = c(1,2), lty = 1,main = "In Sample Fit", ylab = "Values")
        polygon(xxplot,yyplot,col = "lightblue",border = FALSE)
        for (i in which(above_threshold)) {
          rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
               time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
               col = rgb(1, 0, 0, alpha = 0.2), border = NA)
        }
        lines(ts(Yhat_insample),lwd = 2, col = "blue")
        lines(Y_T[-(1:K)],lwd = 1, col = "black")
        legend("bottomleft",fill = c("black","blue","red"),legend = c("Actual Time Series","In Sample Fit","Predicted Regime 2"))
      }else{
        if(Crit == "Entropy"){
          output = Output_list[[which.min(Entropy_list)]]
          P = output$P
          alpha = output$alpha
          zeta_out = output$zeta_t_T
          Y_T = output$Y_T
          Yhat_insample = output$Yhat_insample


          zetaout_ts = ts(zetaout[1,])
          time_vals = time(zetaout_ts)
          plot(zetaout_ts, type = "l", col = "blue", lwd = 2,
               main = "Probability of Regime 2",
               xlab = "Year", ylab = "Probability")
          above_threshold = zetaout_ts > threshold
          for (i in which(above_threshold)) {
            rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
                 time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
                 col = rgb(1, 0, 0, alpha = 0.2), border = NA)
          }
          lines(zetaout_ts, col = "blue", lwd = 2)
          abline(h = threshold, col = "black", lty = 2)


          zetaout_ts_M = ts(t(zetaout))
          plot(zetaout_ts_M,plot.type = "single", col = c(1:N), lwd = 2)
          plot(zetaout_ts_M, col = c(1:N), lwd = 2)





          xplot = time(Y_T[-(1:K)])
          xxplot = c(xplot,rev(xplot))
          yyplot = c(Y_T[-(1:K)],rev(ts(Yhat_insample)))
          plot(cbind(Y_T[-(1:K)],ts(Yhat_insample)),plot.type = "single", col = c("black","blue"),lwd = c(1,2), lty = 1,main = "In Sample Fit", ylab = "Values")
          polygon(xxplot,yyplot,col = "lightblue",border = FALSE)
          for (i in which(above_threshold)) {
            rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
                 time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
                 col = rgb(1, 0, 0, alpha = 0.2), border = NA)
          }
          lines(ts(Yhat_insample),lwd = 2, col = "blue")
          lines(Y_T[-(1:K)],lwd = 1, col = "black")
          legend("bottomleft",fill = c("black","blue","red"),legend = c("Actual Time Series","In Sample Fit","Predicted Regime 2"))
        }
      }
    }
  }
  return(output)
}
