# Two-stage Procrustes rotation

## Usage
```R
#dependencies
source("two_stg_Procrustes.R")
source("promax.R")
source("simplimax.R")
source("Multiple_Starts_nosnow.R")
library(lavaan)
library(dplyr)

#EFA
data("HolzingerSwineford1939")
res <- HolzingerSwineford1939[,c(7:15)] %>%
  fa(nfactors = 3, rotate = "varimax", fm = "ml") 
A <- res$loadings  

#rotation (without multiple starts)
res_twostg <- TwoStgProc(A, 0.5)
res_promax <- promax_tgt(A)
res_simplimax <- simplimax_tgt(A, nrow(A))

#rotation (with multiple starts)
res_twostg_ms <- MULTIPLE_STARTS("TwoStgProc(A, 0.3)", 100)
res_simplimax_ms <- MULTIPLE_STARTS("simplimax_tgt(A, nrow(A))", 100)

#results
round(A %*% res_twostg_ms$solbest$rotmat, 2)
#    [,1]  [,2]  [,3]
#x1  0.04  0.17  0.55
#x2 -0.11  0.04  0.47
#x3  0.02 -0.07  0.63
#x4  0.03  0.78  0.02
#x5  0.03  0.83 -0.06
#x6  0.01  0.75  0.07
#x7  0.69  0.04 -0.15
#x8  0.67 -0.03  0.09
#x9  0.45  0.03  0.33

round(A %*% res_promax$rotmat, 2)
#   [,1] [,2] [,3]
#x1 0.39 0.67 0.26
#x2 0.18 0.49 0.05
#x3 0.15 0.67 0.22
#x4 0.85 0.29 0.19
#x5 0.87 0.22 0.18
#x6 0.83 0.33 0.19
#x7 0.14 0.02 0.69
#x8 0.14 0.25 0.73
#x9 0.24 0.48 0.58

round(A %*% res_simplimax_ms$solbest$rotmat, 2)
#    [,1]  [,2]  [,3]
#x1  0.59 -0.41  0.12
#x2  0.46 -0.19 -0.06
#x3  0.66 -0.18  0.06
#x4  0.06 -0.85  0.19
#x5 -0.02 -0.87  0.19
#x6  0.11 -0.83  0.17
#x7  0.04 -0.15  0.71
#x8  0.28 -0.16  0.69
#x9  0.47 -0.27  0.49
```
