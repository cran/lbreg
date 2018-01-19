#  Negative loglikelihood (nLL)
negll <- function(beta, x, y){
  .y <- y[,1]
  if( ncol(y) == 1 ){  .n <- 1  }else{ .n <- rowSums(y) }
  eta <- x %*% as.matrix(beta)
  llik <- dbinom(.y, size = .n, prob = exp( eta ), log = TRUE ) 
  .value <- -sum(llik)
  return(.value)
} 
