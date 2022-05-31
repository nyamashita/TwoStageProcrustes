#######################################################
#### Regularized regression
#### Naoto Yamashita
#######################################################
library(Rcpp)
library(RcppArmadillo)
sourceCpp("w_update.cpp")

reg_reg <- function(X, #matrix of independent vars.
                    Y, #matrix of dependent var.
                    lambda, #tuning parameter
                    itemax=100, #max of iteration
                    Winit=NULL, #initial value for W
                    cpp=TRUE, #logical; use cpp function (fast) or not
                    eps=1e-7){
  
  p <- ncol(X)
  q <- ncol(Y)
  
  #loss function
  f <- function(W){
    sum((Y - X%*%W)^2) + lambda*sum(abs(W))
  }
  
  if(is.null(Winit)){
    W <- matrix(rnorm(p*q),p,q)
  }else{
    W = Winit
  }
  history <- c()
  
  for(ite in 1:itemax){
    if(cpp){
      W <- W_update(X, Y, lambda)
    }else{
      for(i in 1:nrow(W)){
        for(j in 1:ncol(W)){
          Wols <- solve(t(X)%*%X)%*%t(X)%*%Y
          W[i,j] <- sign(Wols[i,j])*plus(abs(Wols[i,j]) - lambda/(2*sum((X[,i])^2)))
        }
      }
    }
    history[ite] <- f(W)
    if(ite>1){
      if(history[ite-1]-history[ite] < eps){break}
    }
  }

  list(history = history, W = W)
}