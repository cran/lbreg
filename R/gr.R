# Gradient of nLL 
gr <- function(beta, x, y){
  p <- exp( x %*% as.matrix(beta) ) 
  gg <- 0 
  for (i in 1:nrow(x)) {
    if(y[i] == 0) {z <- -exp(log(p[i]) - log(1-p[i]))}else{ z <- 1 }   # avoid division by small number
    gg <- gg + c(x[i,]) * z
  }
  return( -gg )
} 

