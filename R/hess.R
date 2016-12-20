###   Hessian of LL 
hess <- function(beta, x, y){
  p <- exp( x %*% as.matrix(beta) ) 
  z <- p*(y-1) / (1-p)^2
  hh <- 0 
  for (i in 1:nrow(x)) {
    hh <- hh + x[i,] %*% t(x[i,]) * z[i]
  }
  return( hh )
} 

