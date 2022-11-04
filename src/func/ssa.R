require(pracma)
require(tsutils)
ssa = function(X, M= 14){ 
  M = 14
  N = length(X)
  lag = lagmatrix(X, 0:-(M - 1))
  lag = apply(lag, 2, replace_na, 0)
  C = t(lag) %*% lag / N
  e = eigen(C)
  lambda = e$values
  RHO = e$vectors
  PC = lag %*% RHO
  
  RC = zeros(N, M)
  for (m in (1:M)) {
    Z = lagmatrix(PC[,m], 0:-(M - 1))
    Z = apply(Z, 2, replace_na, 0)
    RC[,m] = Z %*% RHO[,m] / M
  }
  return(list("PC" = PC, "RC" = RC))
}

normalize <- function(arr, x, y) {
  m = min(arr)
  range = max(arr) - m
  arr = (arr - m)/range
  
  range2 = y-x
  norm= (arr *range2) +x
}

################################################################################
################################################################################
