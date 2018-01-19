###   Hessian of LL 
hess <- function(beta, x, y){
  p <- exp( x %*% as.matrix(beta) ) 
   .y <- y[,1]
  if( ncol(y) == 1 ){  .n <- 1  }else{ .n <- rowSums(y) }
  z <- -p*(.n-.y) / (1-p)^2
  hh <- 0 
  for (i in 1:nrow(x)) {
    hh <- hh + x[i,] %*% t(x[i,]) * z[i]
  }
  return( hh )
} 

