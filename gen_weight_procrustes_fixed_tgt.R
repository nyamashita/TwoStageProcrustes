### Procrustes target rotation with generalized weighting
### Naoto Yamashita

#fixed target
gen_proc_fixed_tgt <- function(A, # a loading matrix to be rotated
                               T, #target
                               weight_func, # weight function
                               eps = 1e-7
){
  
  p <- nrow(A)
  r <- ncol(A)
  
  #target matrix is provided
  
  #weight matrix
  W <- apply(T,c(1,2),weight_func)
  
  #loss function for GPA
  lf <- function(U){
    sum(((T - A %*% solve(t(U)))*W)^2)
  }
  
  #GPA
  res <- GPA(lf, "oblq", diag(r), itemax = 100, eps = eps)
  U <- res$sol
  
  Ar <- A %*% U
  
  #results
  list(
    rotmat = U,
    rotated = Ar,
    target = T,
    weightmat = W,
    history = res$history
  )
}