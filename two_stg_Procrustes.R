##########################################################
### Two Stage Procrustes Rotation
###  Naoto Yamashita
##########################################################
source("regularized_reg.R")
source("gen_weight_procrustes_fixed_tgt.R")
source("Gradient_Projection_Algorithm_general.R")

weight_func <- function(x){
  1/exp(abs(x))^5
}

TwoStgProc <- function(A, #initial loadings
                       lambda, #tuning parameter
                       stage1 = TRUE, #logical; perform stage1
                       stage2 = TRUE, #logical; perform stage2
                       weight_func = NULL, #function for weighting
                       itemax = 1000, #max of iteration
                       target = T, #target matrix; used only stage1 = FALSE
                       eps = 1e-7){
  
  #stage1: target estimation with lasso penalty
  if(stage1){
    lf <- function(T,S){
      sum((T - A%*%S)^2) + lambda*sum(abs(T))
    }
    history <- c()
    S1 <- GPArotation::Random.Start(ncol(A))
    for(ite in 1:itemax){
      #min. ||AS-T||^2 + \lambda*|T| over T
      T <- reg_reg(X = diag(nrow(A)), 
                   Y = A%*%S1, 
                   lambda = lambda)$W
      #min. ||AS-T||^2 over S
      S1 <- solve(t(GPArotation::targetQ(A, Tmat = S1, Target = T)$Th)) 
      history[ite] <- lf(T,S1)
      if(ite > 1){
        if((history[ite-1] - history[ite]) < eps){
          break
        }
      }
    }
  }else{
    S1 <- diag(rep(1,ncol(A)))
    T <- target
  }
  
  #stage2: procrustes with generalized weighting
  if(stage2){
    res <- gen_proc_fixed_tgt(A%*%S1, T, weight_func)
    S2 <- solve(t(res$rotmat))
    S_out <- S2 %*% S1
    S_out <- S_out %*% diag(diag(t(S_out) %*% S_out)^(-1/2)) #scaling
  }
  
  #return values
  if(stage1 & stage2){ #both stages
    list(
      Tgt = T,
      rotmat_stg1 = S1,
      rotmat_stg2 = S2,
      rotmat = S_out,
      history_stg1 = history,
      history_stg2 = res$history,
      #LOSS_MIN = min(history),
      LOSS_MIN = sum((T - A%*%S_out)^2),
      lambda = lambda
    )
  } else if (!stage1){ #only stage2
    list(
      Tgt = T,
      rotmat_stg1 = S1,
      rotmat_stg2 = S2,
      rotmat = S_out,
      history_stg1 = NULL,
      history_stg2 = res$history,
      LOSS_MIN = min(res$history),
      lambda = lambda
    )
  } else if (!stage2){ #only stage1
    list(
      Tgt = T,
      rotmat_stg1 = S1,
      rotmat_stg2 = NULL,
      rotmat = S1,
      history_stg1 = history,
      history_stg2 = NULL,
      LOSS_MIN = min(history),
      lambda = lambda
    )
  }
}


