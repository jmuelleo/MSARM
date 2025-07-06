
#' @title P_ij_next_function
#' @description P_ij_next_function allows its user to estimate the P matrix of the next iteration as part of the EM algorithm, as described in Hamilton (1994, page 695)
#' @param alpha alpha matrix of the last iteration
#' @param P P matrix of the last iteration
#' @param Y_T Data
#' @param N Number of underlying regimes
#' @param K Lag-order of the AR model
#' @param m Number of time periods the conditional likelihood function conditions on
#'
#' @return Returns a matrix of the dimensions (N x N)
#' @export
#'
#' @examples P_ij_next_function(alpha,P,Y_T,N,K,m)
P_ij_next_function = function(alpha,P,Y_T,N,K,m){
  T = length(Y_T)
  zeta_list = zeta_Yt_function(alpha = alpha,P = P,Y_T = Y_T,N = N,K = K) #optimal inference over the regimes
  zeta_t_t = zeta_list[[1]]
  zeta_t1_t = zeta_list[[2]]
  zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K) #smoothed inference over the regimes

  P_i_j_matrix_t_function = function(alpha,P,t,Y_T,N,K){
    #Creates the matrix of probabilites for going from i to j for the time points t-1 to t
    P_i_j_t = matrix(5,ncol = N,nrow = N)
    for(i in 1:N){ #Order of filling the matrix: 1,1; 1,2; 2,1; 2,2
      for(j in 1:N){
        P_i_j_t[j,i] = zeta_YT[j,t]*((P[j,i]*zeta_t_t[i,t-1])/zeta_t1_t[j,t-1])
      }
    }
    return(P_i_j_t)
  }

  Cumsum_matrix = matrix(0,ncol = N,nrow = N)
  for(t in (m+1+K):T){ #Here the conditional log-likelihood is optimised, m is the number of time periods one conditions on
    matrix_add = P_i_j_matrix_t_function(alpha,P,t,Y_T,N,K)
    Cumsum_matrix = Cumsum_matrix + matrix_add
  }

  sum_ij_gemeinsam = Cumsum_matrix

  Pnext = matrix(5,ncol = N,nrow = N)
  for(i in 1:N){
    sum_i = sum(zeta_YT[i,(K+1+m-1):(T-1)]) #-1 on both sides, as s_{t-1} = i is of interest
    Pnext[,i] = sum_ij_gemeinsam[,i]/sum_i #i to j, therefore ith column divided by sum_i
  }

  #Is done to reduce rounding errors, normally the rounding errors are around the size 0.002
  Pnext_adj = matrix(5,ncol = N,nrow = N)
  for(i in 1:N){
    Pnext_adj[,i] = Pnext[,i]/apply(Pnext,2,sum)[i]
  }


  return(Pnext_adj)

}
