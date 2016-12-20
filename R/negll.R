#  Negative loglikelihood (nLL)
negll <- function(beta, x, y){
  eta <- x %*% as.matrix(beta)
  llik <- dbinom(y, size = 1, prob = exp( eta ), log = TRUE ) 
  .value <- -sum(llik)
  return(.value)
} 
