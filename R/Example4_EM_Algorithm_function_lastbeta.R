
#' @title Example4_EM_Algorithm_function_lastbeta
#' @description Example4_EM_Algorithm_function_lastbeta allows its user to run the EM algorithm for estimating the parameters of a Markov-Switching AR model of the form of Example 4, see Müller (2025, page 19-21)
#' @param Y_T Data
#' @param K Lag-order of the AR model
#' @param N Number of underlying regimes
#' @param m Number of time periods the conditional log-likelihood conditions on
#' @param threshold Threshold for assigning the regimes
#' @param max Number of iterations of the EM algorithm
#' @param R Number of random starting points of the EM algorithm
#' @param Crit Performance metric used for choosing between the optima
#' @param all.plot Whether the plots from all optimisation attempts should be printed
#' @param Switcher Switching vector indicating for all parameters whether they switch or not
#'
#' @return Returns a list object containing the results from the EM algorithm estimation of the Markov-Switching AR model
#' @export
#'
#' @examples Example4_EM_Algorithm_function_lastbeta(Y_T,K,N,m,threshold = 0.5,max = 500,R = 15, Crit = "LV", all.plot = FALSE, Switcher)
Example4_EM_Algorithm_function_lastbeta = function(Y_T,K,N,m,threshold = 0.5,max = 500,R = 15, Crit = "LV", all.plot = FALSE, Switcher){

  Output_list = as.list(1:R)

  LV_list = as.list(1:R)
  RSS_list = as.list(1:R)
  RCM_list = as.list(1:R)
  Entropy_list = as.list(1:R)
  for(r in 1:R){


    zetahat_1_0 = zeta_1_0_function(Y_T,N,K) #random starting point of the EM algorithm
    alpha_0 = alpha_0_function(Y_T,N,K) #random starting point of the EM algorithm
    P_0 = P_0_function(Y_T,N,K) #random starting point of the EM algorithm
    T = length(Y_T)
    counter = 1
    P_Sammler = as.list(1:(max+1))
    alpha_Sammler = as.list(1:(max+1))
    P_Sammler[[1]] = P_0
    alpha_Sammler[[1]] = alpha_0
    #run the EM algorithm for max iterations
    while(counter <= max){

      P_Sammler[[counter+1]] = P_ij_next_function(alpha = alpha_Sammler[[counter]],P = P_Sammler[[counter]],Y_T = Y_T,N = N,K = K,m) #Computing the P matrix of the next iteration
      alpha_Sammler[[counter+1]] = Example4_alpha_next_function_lastbeta(Y_T,alpha_Sammler[[counter]],P_Sammler[[counter]],Switcher,N,K,m) #Computing the alpha matrix of the next iteration
      Log_value_l = Log_Likelihood_function(alpha_Sammler[[ifelse(counter-1>=1,counter-1,1)]],P_Sammler[[ifelse(counter-1>=1,counter-1,  1)]],Y_T,N,K) #Value of the log-likelihood from the last iteration
      Log_value_n = Log_Likelihood_function(alpha_Sammler[[counter]],P_Sammler[[counter]],Y_T,N,K) #value of the log-likelihood from the current iteration

      Delta = abs(Log_value_n - Log_value_l) #Change in the value of the log-likelihood

      # if(Delta <= 0.00001){  #Several experiments showed that simply choosing a fixed number of iterations leads to a better result than setting a minimum change in the value of the log-likelihood, still this approach might be reimplemented in future versions of MSARM
      #   max = counter - 1
      #   break
      #}

      counter = counter + 1
    }
    Log_value = Log_Likelihood_function(alpha_Sammler[[counter]],P_Sammler[[counter]],Y_T,N,K) #value of the conditional likelihood based on the final estimates



    zetaout = zeta_YT_function(alpha_Sammler[[max+1]],P_Sammler[[max+1]],Y_T,N,K)[,-(1:(K))] #smoothed inference over the regimes
    zetaout_ts = ts(zetaout[1,]) #1 oder 2 wählen je nachdem, wie das Modell die Regime zuordnet

    time_vals = time(zetaout_ts)
    above_threshold = zetaout_ts > threshold
    #plot is printed if all.plot = TRUE
    if(all.plot == TRUE){
    plot(zetaout_ts, type = "l", col = "blue", lwd = 2,
         main = "Probability of Regime 2",
         xlab = "Year", ylab = "Probability")
    for (i in which(above_threshold)) {
      rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
           time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
           col = rgb(1, 0, 0, alpha = 0.2), border = NA)
    }
    lines(zetaout_ts, col = "blue", lwd = 2)
    abline(h = threshold, col = "black", lty = 2)
}
    years = floor(time(zetaout_ts))

    zetaout_ts_M = ts(t(zetaout))
    #plot is printed if all.plot = TRUE
if(all.plot == TRUE){
    plot(zetaout_ts_M,plot.type = "single", col = c(1:N), lwd = 2)
    plot(zetaout_ts_M, col = c(1:N), lwd = 2)
}
    Predicted_Regime = apply(zetaout_ts_M,1,function(row){which(row == max(row))}) #Predicted Regimes

    Yhat_insample = rep(0,T-K)
    for(t in 1:(T-K)){
      for(j in 1:N){

        if(j == Predicted_Regime[t]){

          Yhat_insample[t] = c(c(1,Y_T[(t+K-1):(t)])%*%alpha_Sammler[[max+1]][j,-(K+2)]) #Computing in-sample fit

        }

      }
    }
    #plot is printed if all.plot = TRUE
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
    legend("bottomleft",fill = c("black","blue","red"),legend = c("Actual Time Series","In Sample Fit","Predicted Regime"))
}
    #Computing RSS as performance metric measuring in-sample fit
    RSS = sum((Y_T[-(1:K)] - ts(Yhat_insample))^2)
    residuals = Y_T[-(1:K)] - Yhat_insample

    #Computing the RCM value
    RCM = 100*(N^2)*mean(apply(zetaout,2,prod))


    Persistency = prod(diag(P_Sammler[[max+1]])) #Earlier versions implemented this as additional performance metric, but did not lead to better results, might be reimplemented in an adjusted form in future MSARM updates

    #Computes the Entropy value
    Entropy = 100*(N^2)*(-1)*mean(apply(zetaout,2,function(col){
      phat = col[1]
      return(phat*log(phat))
    }))

    Entropy = Entropy

    #Saving the results
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
  #Create the output and plot for the optimum selected by LV
  if(Crit == "LV"){
    output = Output_list[[which.max(LV_list)]]
    P = output$P
    alpha = output$alpha
    zetaout = output$zeta_t_T
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


    regime_colors <- sapply(1:(N - 1), function(l) {
      rgb(1 / l, ifelse(l > 1, 1 / (l - 1), 0), ifelse(l > 2, 1 / (l - 2), 0), alpha = 0.2)
    })

    zetout_length = length(zetaout_ts)
    above_threshold_matrix = matrix(0,ncol = zetout_length, nrow = N-1 )
    for(l in 1:(N-1)){
      above_threshold_matrix[l,] = ts(zetaout[l,]) > threshold


      for (i in which(above_threshold_matrix[l,] == 1)) {
        rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
             time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
             col = rgb(1/l, ifelse(l > 1, 1/(l-1), 0), ifelse(l > 2, 1/(l-2), 0), alpha = 0.2), border = NA)
      }
    }



    lines(ts(Yhat_insample),lwd = 2, col = "blue")
    lines(Y_T[-(1:K)],lwd = 1, col = "black")

    legend_labels <- c("Actual Time Series", "In Sample Fit", paste0("Regime ", 1:(N - 1)))
    legend_colors <- c("black", "blue", regime_colors)
    legend("bottomleft", fill = legend_colors, legend = legend_labels)
  }else{
    #Create the output and plot for the optimum selected by RSS
    if(Crit == "RSS"){
      output = Output_list[[which.min(RSS_list)]]
      P = output$P
      alpha = output$alpha
      zetaout = output$zeta_t_T
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


      regime_colors <- sapply(1:(N - 1), function(l) {
        rgb(1 / l, ifelse(l > 1, 1 / (l - 1), 0), ifelse(l > 2, 1 / (l - 2), 0), alpha = 0.2)
      })

      zetout_length = length(zetaout_ts)
      above_threshold_matrix = matrix(0,ncol = zetout_length, nrow = N-1 )
      for(l in 1:(N-1)){
        above_threshold_matrix[l,] = ts(zetaout[l,]) > threshold


        for (i in which(above_threshold_matrix[l,] == 1)) {
          rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
               time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
               col = rgb(1/l, ifelse(l > 1, 1/(l-1), 0), ifelse(l > 2, 1/(l-2), 0), alpha = 0.2), border = NA)
        }
      }



      lines(ts(Yhat_insample),lwd = 2, col = "blue")
      lines(Y_T[-(1:K)],lwd = 1, col = "black")

      legend_labels <- c("Actual Time Series", "In Sample Fit", paste0("Regime ", 1:(N - 1)))
      legend_colors <- c("black", "blue", regime_colors)
      legend("bottomleft", fill = legend_colors, legend = legend_labels)
    }else{
      #Create the output and plot for the optimum selected by RCM
      if(Crit == "RCM"){
        output = Output_list[[which.min(RCM_list)]]
        P = output$P
        alpha = output$alpha
        zetaout = output$zeta_t_T
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


        regime_colors <- sapply(1:(N - 1), function(l) {
          rgb(1 / l, ifelse(l > 1, 1 / (l - 1), 0), ifelse(l > 2, 1 / (l - 2), 0), alpha = 0.2)
        })

        zetout_length = length(zetaout_ts)
        above_threshold_matrix = matrix(0,ncol = zetout_length, nrow = N-1 )
        for(l in 1:(N-1)){
          above_threshold_matrix[l,] = ts(zetaout[l,]) > threshold


          for (i in which(above_threshold_matrix[l,] == 1)) {
            rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
                 time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
                 col = rgb(1/l, ifelse(l > 1, 1/(l-1), 0), ifelse(l > 2, 1/(l-2), 0), alpha = 0.2), border = NA)
          }
        }



        lines(ts(Yhat_insample),lwd = 2, col = "blue")
        lines(Y_T[-(1:K)],lwd = 1, col = "black")

        legend_labels <- c("Actual Time Series", "In Sample Fit", paste0("Regime ", 1:(N - 1)))
        legend_colors <- c("black", "blue", regime_colors)
        legend("bottomleft", fill = legend_colors, legend = legend_labels)
      }else{
        #Create the output and plot for the optimum selected by Entropy
        if(Crit == "Entropy"){
          output = Output_list[[which.min(Entropy_list)]]
          P = output$P
          alpha = output$alpha
          zetaout = output$zeta_t_T
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


          regime_colors <- sapply(1:(N - 1), function(l) {
            rgb(1 / l, ifelse(l > 1, 1 / (l - 1), 0), ifelse(l > 2, 1 / (l - 2), 0), alpha = 0.2)
          })

          zetout_length = length(zetaout_ts)
          above_threshold_matrix = matrix(0,ncol = zetout_length, nrow = N-1 )
          for(l in 1:(N-1)){
            above_threshold_matrix[l,] = ts(zetaout[l,]) > threshold


            for (i in which(above_threshold_matrix[l,] == 1)) {
              rect(time_vals[i] - 0.5 / frequency(zetaout_ts), par("usr")[3],
                   time_vals[i] + 0.5 / frequency(zetaout_ts), par("usr")[4],
                   col = rgb(1/l, ifelse(l > 1, 1/(l-1), 0), ifelse(l > 2, 1/(l-2), 0), alpha = 0.2), border = NA)
            }
          }



          lines(ts(Yhat_insample),lwd = 2, col = "blue")
          lines(Y_T[-(1:K)],lwd = 1, col = "black")

          legend_labels <- c("Actual Time Series", "In Sample Fit", paste0("Regime ", 1:(N - 1)))
          legend_colors <- c("black", "blue", regime_colors)
          legend("bottomleft", fill = legend_colors, legend = legend_labels)
        }
      }
    }
  }
  return(output)
}
