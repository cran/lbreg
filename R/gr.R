# Gradient of nLL 
gr <- function(beta, x, y){
  p <- exp( x %*% as.matrix(beta) ) 
  .y <- y[,1]
  gg <- 0 
  
  if( ncol(y) == 1 ){  
  
    for (i in 1:nrow(x)) {
      if(.y[i] == 0) {z <- -exp(log(p[i]) - log(1-p[i]))}else{ z <- 1 }   # avoid division by small number
      gg <- gg + c(x[i,]) * z
    }
  }else{ 
	.n <- rowSums(y) 
	for (i in 1:nrow(x)) {
	  z <- (.y[i] - .n[i]*p[i])/(1-p[i])  
	  gg <- gg + c(x[i,]) * z
    }
  }
  return( -gg )
} 

