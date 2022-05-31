##########################################################
### Gradient Projection Algorithm
###  Naoto Yamashita
##########################################################

GPA <- function(LF, #function to be minimized
                method, #orthgonal or oblique
                Qini, #initial value for parameter matrix
                itemax=1000, #max of iteration
                eps=0.0001,
                verbose=FALSE){
  
  
  #function for numarical derivative
  #evaluate numerical derivative at current Q
  Nudev <- function(Q){
    Jmat <- function(row,col){
      J <- matrix(0,nrow(Q),ncol(Q))
      J[row,col] <- 1
      J
    }
    G <- matrix(0,nrow(Q),ncol(Q))
    for(row in 1:nrow(Q)){
      for(col in 1:ncol(Q)){
        G[row,col] <- (LF(Q+eps*Jmat(row,col)) - LF(Q-eps*Jmat(row,col)))/(2*eps)
      }
    }
    G
  }
  
  ##LF is minimized over rotation matrix Q
  ##LF must given by a MINIMIZATION criterion
  
  #update function for oblique
  Qup_oblq <- function(Q,Gra,alpha){
    mat <- Q - alpha*Gra
    mat %*% solve(sqrt(diag(diag(t(mat)%*%mat))))
  }
  
  #update function for otrhogonal
  Qup_orth <- function(Q,Gra,alpha){
    mat <- Q - alpha*Gra
    res <- svd(mat)
    res$u %*% t(res$v)
  }
  
  #initial values
  Q <- Qini
  alphaini <- 1
  alpha <- alphaini
  history <- c()
  
  #iteration starts
  for(ite in 1:itemax){
    #compute gradient matrix
    Gra <- Nudev(Q)
    repeat{
      #update Q
      if(method == "orth"){
        Q2 <- Qup_orth(Q,Gra,alpha)
      }
      if(method == "oblq"){
        Q2 <- Qup_oblq(Q,Gra,alpha)
      }
      #if current alpha decreases LF()...
      if(LF(Q) > LF(Q2)){break}
      
      alpha <- alpha/2
      if(alpha < 0.0001){break}
    } 
    Q <- Q2
    history[ite] <- LF(Q)
    #initialize alpha
    alpha <- alphaini
    #stopping rule
    if(ite > 1){
      if((history[ite - 1] - history[ite]) < 0.0000001){break} 
    }
    if(ite == itemax & verbose){
      print("Algorithm did not converged. Step size of numerical derivative(eps) should be reconsidered.")
    }     
  }
  list(history=history,sol=Q, LOSS_MIN=min(history))
}