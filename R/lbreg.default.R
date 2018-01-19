lbreg.default <- function(x, y, start.beta, tol=.9999, delta=1, ...)
{
  y <- as.matrix(y)	
  if( ncol(y) == 1 ){  .y <- y  }else{ .y <- y[,1] }
  if( missing(start.beta) ){
	b0 <- coef( glm.fit(x, .y, family = poisson(link = "log") ) ) 
	mX <- as.matrix(-x[,-1])
	b0[1] <-  min( mX %*% b0[-1] ) - delta
  }else{
	b0 <- start.beta
  }

  fit <- lbreg.fit(x=as.matrix(x), y=y, start.beta=b0, tol=tol)  

  fit$call <- match.call()
  return(fit)  

}


