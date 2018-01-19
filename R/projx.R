
projx <- function(x){
  u = svd(x)$u
  return(u%*%t(u))
}

