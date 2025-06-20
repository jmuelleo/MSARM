#' Title
#'
#' @param alpha alpha vector (intercept, coefficients and error variance) from the last iteration
#' @param P Transpose of the transition matrix Pi
#' @param Y_T Data
#' @param N Number of underlying Regimes
#' @param K Lag order of the AR-Model
#'
#' @return Computes the P matrix for the next iteration, based on the last iteration
#' @export
#'
#' @examples P_ij_next_function(alpha,P,Y_T,N,K)
P_ij_next_function = function(alpha,P,Y_T,N,K,m){
  T = length(Y_T)
  zeta_list = zeta_Yt_function(alpha = alpha,P = P,Y_T = Y_T,N = N,K = K)
  zeta_t_t = zeta_list[[1]]
  zeta_t1_t = zeta_list[[2]]
  zeta_YT = zeta_YT_function(alpha,P,Y_T,N,K)

  P_i_j_matrix_t_function = function(alpha,P,t,Y_T,N,K){
    #bestimmt die Matrix der gemeinsamen Wahrscheinlichkeiten i nach j zum Zeitpunkt t-1 zu t
    P_i_j_t = matrix(5,ncol = N,nrow = N)
    for(i in 1:N){ #Hier ist die Schleifenreihenfolge für die Befüllung der Matrix wichtig 1,1; 1,2; 2,1; 2,2
      for(j in 1:N){
        P_i_j_t[j,i] = zeta_YT[j,t]*((P[j,i]*zeta_t_t[i,t-1])/zeta_t1_t[j,t-1])
      }
    }
    return(P_i_j_t)
  }

  Cumsum_matrix = matrix(0,ncol = N,nrow = N)
  for(t in (m+1+K):T){ #Hier wird die bedingte Log-Likelihood optimiert, das ist m, die Menge an Bedingungszeitpunkten
    matrix_add = P_i_j_matrix_t_function(alpha,P,t,Y_T,N,K)
    Cumsum_matrix = Cumsum_matrix + matrix_add
  }

  sum_ij_gemeinsam = Cumsum_matrix

  Pnext = matrix(5,ncol = N,nrow = N)
  for(i in 1:N){
    sum_i = sum(zeta_YT[i,(K+1+m-1):(T-1)]) #-1 auf beiden Seiten, da s_{t-1} = i betrachtet wird
    Pnext[,i] = sum_ij_gemeinsam[,i]/sum_i #i nach j, also muss die ite Spalte durch sum_i geteilt werden
  }

  Pnext_adj = matrix(5,ncol = N,nrow = N)
  for(i in 1:N){
    Pnext_adj[,i] = Pnext[,i]/apply(Pnext,2,sum)[i]
  }
  #Einfach um kleine Rundungsfehler rauszunehmen, so dass sich die Spalten genau auf 1 summieren (Größenordnung 0.002 ca.)

  return(Pnext_adj)

}
