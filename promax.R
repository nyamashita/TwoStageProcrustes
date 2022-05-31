########################################################
#### Promax Rotation
####  Naoto Yamashita
########################################################

library(GPArotation)

promax_tgt <- function(
  A, #loading matrix
  alpha = 4, #tuning parameter
  normalize = TRUE #Kaiser's normalization
){
  if(normalize){
    A <- diag(diag(A%*%t(A))^(-1/2)) %*% A
  }
  
 #varimax
  varimax_A <- A %*% varimax(A)$rotmat
  Tgt <- abs(A) * varimax_A^{alpha-1}
  
  #target rotation
  res <- targetQ(A, Target = Tgt)
  
  #return
  list(
    rotmat = res$Th,
    loadings = A %*% t(solve(res$Th)),
    Tgt = Tgt,
    convergence = res$convergence
  )
}