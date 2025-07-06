
#' @title Example3_EM_Algorithm_function_lastbeta
#' @description Example3_EM_Algorithm_function_lastbeta allows its user to apply the EM algorithm for estimating Markov-Switching models of the form of Example 3 (All Switching), as described in Müller (2025, page 18-19)
#' @param Y_T Data
#' @param K Lag-order of the AR model
#' @param N Number of underlying regimes
#' @param m Number of time periods the conditional likelihood conditions on
#' @param threshold #threshold for assigning the predicted regimes
#' @param max #number of iterations for the EM algorithm
#' @param R #number of random starting points for the EM algorithm
#' @param Crit #performance metric for choosing one of the R optima
#' @param all.plot #whether plots from all R optimisation attempts should be printed
#' @param Switcher #Switching vector indicating which parameters switch
#'
#' @return returns a list object containg the results from the EM algorithm estimation of the Markov-Switching AR model
#' @export
#'
#' @examples Example3_EM_Algorithm_function_lastbeta(Y_T,K,N,m,threshold = 0.5,max = 500,R = 15, Crit = "LV", all.plot = FALSE, Switcher)
Example3_EM_Algorithm_function_lastbeta = function(Y_T,K,N,m,threshold = 0.5,max = 500,R = 15, Crit = "LV", all.plot = FALSE, Switcher){
  T = length(Y_T)
  Output_list = as.list(1:R)

  LV_list = as.list(1:R)
  RSS_list = as.list(1:R)
  RCM_list = as.list(1:R)
  Entropy_list = as.list(1:R)
  for(r in 1:R){


    zetahat_1_0 = zeta_1_0_function(Y_T,N,K) #random starting points
    alpha_0 = alpha_0_function(Y_T,N,K) #random starting points
    P_0 = P_0_function(Y_T,N,K) #random starting points

    counter = 1 #begin EM algorithm
    P_Sammler = as.list(1:(max+1))
    alpha_Sammler = as.list(1:(max+1))
    P_Sammler[[1]] = P_0
    alpha_Sammler[[1]] = alpha_0
    while(counter <= max){ #run the EM algorithm until iteration max is reached

      P_Sammler[[counter+1]] = P_ij_next_function(alpha_Sammler[[counter]],P_Sammler[[counter]],Y_T,N,K,m) #compute next P matrix
      alpha_Sammler[[counter+1]] = Example3_alpha_next_function(alpha_Sammler[[counter]],P_Sammler[[counter]],Y_T,N,K,m) #compute next alpha matrix
      Log_value_l = Log_Likelihood_function(alpha_Sammler[[ifelse(counter-1>=1,counter-1,1)]],P_Sammler[[ifelse(counter-1>=1,counter-1,  1)]],Y_T,N,K) #last log likelihood value
      Log_value_n = Log_Likelihood_function(alpha_Sammler[[counter]],P_Sammler[[counter]],Y_T,N,K) #current log likelihood value

      Delta = abs(Log_value_n - Log_value_l) #change in log likelihood value

      #  if(Delta <= 0.00001){ #several experiments have revealed that just doing a pre set number of iterations leads to better results than setting a minimum change, but later package versions might utilize this
      #    max = counter - 1
      #    break
      # }

      counter = counter + 1
    }
    Log_value = Log_Likelihood_function(alpha_Sammler[[counter]],P_Sammler[[counter]],Y_T,N,K) #value of the log likelihood





    zetaout = zeta_YT_function(alpha_Sammler[[max+1]],P_Sammler[[max+1]],Y_T,N,K)[,-(1:(K))] #smoothed inference based on the parameter results
    zetaout_ts = ts(zetaout[1,])

    #plotting the results if all.plot = TRUE
    time_vals = time(zetaout_ts)
    if(all.plot == TRUE){
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

#plotting results if all.plot = TRUE
    zetaout_ts_M = ts(t(zetaout))
    if(all.plot == TRUE){
    plot(zetaout_ts_M,plot.type = "single", col = c(1:N), lwd = 2)
    plot(zetaout_ts_M, col = c(1:N), lwd = 2)}

    #computing the predicted regimes
    Predicted_Regime = apply(zetaout_ts_M,1,function(row){which(row == max(row))})

    #creating the in-sample fit
    Yhat_insample = rep(0,T-K)
    for(t in 1:(T-K)){
      for(j in 1:N){

        if(j == Predicted_Regime[t]){

          Yhat_insample[t] = c(c(1,Y_T[(t+K-1):(t)])%*%alpha_Sammler[[max+1]][j,-(K+2)])

        }

      }
    }

    #plot of the in-sample results if all.plot = TRUE
    xplot = time(Y_T[-(1:K)])
    xxplot = c(xplot,rev(xplot))
    yyplot = c(Y_T[-(1:K)],rev(ts(Yhat_insample)))
    if(all.plot == TRUE){
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

    #computing RSS as in-sample fit metric
    RSS = sum((Y_T[-(1:K)] - ts(Yhat_insample))^2)
    residuals = Y_T[-(1:K)] - Yhat_insample

    #computing RCM
    RCM = 100*(N^2)*mean(apply(zetaout,2,prod))


    Persistency = prod(diag(P_Sammler[[max+1]])) #earlier attempts utilized this performance metric, later versions of MSARM might reimplement such an approach

    #computing the Entropy value
    Entropy = 100*(N^2)*(-1)*mean(apply(zetaout,2,function(col){
      phat = col[1]
      return(phat*log(phat))
    }))

    Entropy = Entropy

    #saving the results
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




#Choosing the optima based on LV and printing the results and plots
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












#Choosing the optima based on RSS and printing the results and plots
  }else{
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


#Choosing the optima based on RCM and printing the results and plots
    }else{
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


#Choosing the optima based on Entropy and printing the results and plots
      }else{
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
