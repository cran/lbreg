predict.lbreg <- function(object, newdata, ...)
{
  #  prediction
   if(missing(newdata)){
     return( fitted(object) )
   }else{
	  mf <- delete.response(terms(object))
      X <- model.frame(mf, data=object$data)
      Z <- model.frame(mf, data=newdata) 
      Xz <- rbind(X, Z) 
      X <- model.matrix(terms(mf), data=X)
      Xz <- model.matrix(terms(mf), data=Xz)
        
      y <- as.matrix( model.response(model.frame(object)) )
   
      np <- nrow(newdata)
      n <- nrow(X) 
      
      out.z <- constrOptim(theta = object$start.beta, f = negll, grad = gr, method="BFGS", 
            ui = -Xz, ci = rep(0,nrow(Xz)),
            x = X, y = y)
      beta.pred <- out.z$par
      yhat <- as.vector(exp(Xz %*% beta.pred))
      ypred <- yhat[-(1:n)]                       # not efficient 
  
      J <- -hess(beta=beta.pred, x=X, y=y)   # Obs Fisher Info
    
      iJ <- MASS::ginv( J )
   # active constraints
   a <- which( yhat >= object$tol )
   Afull <- NULL 
   if( any(a) ){     # check for active constraints and proceed with correction of VCOV 
	Afull = Xz[a,,drop=FALSE]
    A = Afull
    A = unique(A) 
    r = NROW(A)
    p = NCOL(A)
    
    if(r < p){
    Wtil = diag(p) - projx(t(A))
    W = Wtil[1:p,1:(p-r)]
    
    if( max( A%*%W ) > sqrt( 1e-6 ) ) { warning("max AW is > 1e-6") }
    
    meat <- MASS::ginv(t(W) %*% J %*% W) 
    V <- W %*% meat %*% t(W)
  }else{
	V <- iJ - iJ %*% t(A) %*% MASS::ginv(A %*% iJ %*% t(A)) %*% A %*% iJ
   }
    #nnew <- nrow(Xz) - nrow(X)
    #rownames(Afull) <- c(1:nrow(X), paste('new', 1:nnew, sep=''))[a]
    #colnames(Afull) <- colnames(X)
        # no active constraints
    }else{   
		V <- iJ
   }

    # SE of prediction
    
    se.pred <- double(np)
    
    for(i in 1:np){
	  z <- c( Xz[n+i,] )
      gp <- z*ypred[i]    # derivative of exp( linpred )
      dim(gp) <- c(length(gp),1)
      se.pred[i] <- sqrt( t(gp) %*% V %*% gp ) # delta method
    } 
    
        
   return(list(
       ypred = ypred,
       se.pred=se.pred, 
       coef.pred = beta.pred, 
       convergence = out.z$convergence,
       Active = Afull, 
       tol = object$tol
       )) 
   }
}
