predict.lbreg <- function(object, newdata, ...)
{
  #  prediction
   if(is.null(newdata)){
     return( fitted(object) )
   }else{
      XX <- model.frame(formula=object$formula, data=object$data)
      X <- model.matrix(attr(XX, "terms"), data=XX)
      
      np <- nrow(newdata)
      n <- nrow(object$data)
      
      novo <- cbind(0.5, newdata)    # couldn't use NA because model.frame below omits rows with NA; use 0.5 !!!!
      colnames(novo)[1] <- as.character( object$formula[[2]] )   
      NewData <- rbind(XX,novo)
      
      mf <- model.frame(formula=object$formula, data=NewData)
      Xz <- model.matrix(attr(mf, "terms"), data=mf)
      
      y <- model.response(mf)[1:n]
 
      out.z <- constrOptim(theta = object$start.beta, f = negll, grad = gr, method="BFGS", 
            ui = -Xz, ci = rep(0,nrow(Xz)),
            x = X, y = y)
      beta.pred <- out.z$par
      yhat <- as.vector(exp(Xz %*% beta.pred))
      ypred <- yhat[-(1:n)]                       # not efficient 
  
   J <- -hess(beta=beta.pred, x=X, y=y)   # Obs Fisher Info
    
   # active constraints
   a <- yhat >= object$tol
   A <- active <- NULL 
   if( any(a) ){     # check for active constraints and proceed with correction of VCOV 
	Afull = Xz[a,]
    A = Afull
    if(sum(a) == 1){
     r = 1
     p = ncol(Xz)
     dim(A) = c(r,p)
    }else{
     A = unique(A) 
     r = nrow(A)
     p = ncol(Xz)
     dim(A) = c(r,p)
    }
    Wtil = diag(p) - proj(t(A))
    W = Wtil[1:p,1:(p-r)]
    if( max( A%*%W ) > sqrt( 1e-6 ) ) { warning("max AW is > 1e-6") }
    meat <- try( MASS::ginv(t(W) %*% J %*% W) )
      if(class(meat) != "try-error"){
         V <- W %*% meat %*% t(W)
      }else{
         stop("could not invert (W'JW) --- 'vcov' ignores active constraints") 
    }
 
  
    active <- Afull
    nnew <- nrow(Xz) - nrow(X)
    rownames(active) <- c(1:nrow(X), paste('new', 1:nnew, sep=''))[a]
    colnames(active) <- colnames(X)
    }else{   # no active constraints
		V <- MASS::ginv(J)
	}

    # SE of prediction
    
    se.pred <- double(np)
    
    for(i in 1:np){
	  z <- c( Xz[n+i,] )
      gp <- z*ypred    # derivative of exp( linpred )
      dim(gp) <- c(length(gp),1)
      se.pred[i] <- sqrt( t(gp) %*% V %*% gp ) # delta method
    } 
    
        
   return(list(
       ypred = ypred,
       se.pred=se.pred, 
       coef.pred = beta.pred, 
       convergence = out.z$convergence,
       Active.pred = active, 
       tol = object$tol
       )) 
   }
}
