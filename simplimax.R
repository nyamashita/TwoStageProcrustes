########################################################
#### Simplimax Rotation
####  Naoto Yamashita
########################################################

library(GPArotation)

simplimax_tgt <- function(
  A, #loading matrix
  k = (ncol(A)-1)*nrow(A), #tuning parameter; number of zeros in T
  itemax = 100,
  eps = 1e-7
){
  #loss function
  lf <- function(U, T){
    sum((A%*%t(solve(U)) - T)^2)
  }
  
  #initialization
  U <- Random.Start(ncol(A))
  
  history <- c()
  for(ite in 1:itemax){
    #T-step
    Ttmp <- A%*%t(solve(U))
    T_base <- abs(Ttmp) > sort(abs(Ttmp))[k]
    T <- T_base*Ttmp
    
    #U-step
    U <- targetQ(A, Target = T, maxit = 10)$Th
    
    #check convergence
    history[ite] <- lf(U, T)
    if(ite > 1){
      if(history[ite-1] - history[ite] < eps){
        break
      }
    }
  }
  
  #return
  list(
    rotmat = U,
    loadings = A %*% t(U),
    Tgt = T,
    history = history,
    LOSS_MIN = min(history)
  )
}
