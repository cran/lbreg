lbreg.default <- function(x, y, start.beta, tol=.9999, delta=1, ...)
{
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  if( missing(start.beta) ){
	b0 <- coef( glm.fit(x, y, family = poisson(link = "log") ) ) 
	mX <- as.matrix(-x[,-1])
	b0[1] <-  min( mX %*% b0[-1] ) - delta
  }else{
	b0 <- start.beta
  }

  est <- lbregEst(x, y, start.beta=b0, tol)  

  est$residuals <- y - est$fitted.values

  est$call <- match.call()

  class(est) <- "lbreg"
  return(est)  
}


