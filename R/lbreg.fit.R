lbreg.fit <- function(x, y, start.beta, tol, ...)
{
    X <- x
    #y <- as.matrix(y)  
    out <- constrOptim(theta = start.beta, f = negll, grad = gr, method="BFGS", hessian=TRUE,
            ui = -X, ci = rep(0,nrow(X)),
            x = X, y = y, control = list(maxit=200))

    beta <- out$par  # MLE
    
    ### Uncorrected Var(mle)
    #vcov0 <- MASS::ginv( out$hess )   # using numerical hessian from optimization at MLE
    J <- -hess(beta=beta, x=X, y=y)   # Obs Fisher Info
    vcov0 <- MASS::ginv( J )          # using our function hess at MLE  
    colnames(vcov0) <- rownames(vcov0) <- names(beta) <- colnames(X)  
    ###  end

 
    yhat <- as.vector(exp(X %*% beta))
    #out$fitted.values <- yhat
    #out$residuals <- y - yhat
    
   # active constraints
   a <- yhat >= tol
   A <- active <- NULL 
   if( any(a) ){     # check for active constraints and proceed with correction of VCOV 
    Afull = X[a,]
    A = Afull
    if(sum(a) == 1){
     r = 1
     p = ncol(X)
     dim(A) = c(r,p)
    }else{
     A = unique(A)   # Afull may have repeated rows
     r = nrow(A)
     p = ncol(X)
     dim(A) = c(r,p)
    }
    Wtil = diag(p) - projx(t(A))
    W = Wtil[1:p,1:(p-r)]
    if( max( A%*%W ) > sqrt( 1e-6 ) ) { warning("max AW is > 1e-6") }
    meat = try( MASS::ginv(t(W) %*% J %*% W) )
      if(class(meat) != "try-error"){
         new.vcov <- W %*% meat %*% t(W)
         out$vcov <- new.vcov
      }else{
         warning("could not invert (W'JW) --- reported 'vcov' is ignoring active constraints")
    }
   active <- Afull
   if( any(a) ){  dim(active) <- dim(A)  }
   rownames(active) <- (1:nrow(X))[a]
   colnames(active) <- colnames(X)
   }else{  # no active constraints
	   out$vcov <- vcov0
   }   

   SE <- sqrt( diag( out$vcov ) )
   names(SE) <- names(beta)
   
   rownames(out$vcov) <- colnames(out$vcov) <- names(beta)
  
   # deviance residuals 
   y1 <- y[,1]
   if( ncol(y) == 1 ){  
	     n <- 1
	     i <- y1 == 0 
         j <- y1 == 1 
         Y <- y1
         d <- y1*log(Y/yhat) + (1-Y)*log((1-Y)/(1-yhat)) 
	     di <- -log(1-yhat[i])
	     dj <- -log(yhat[j])
	     d[i] <- di
	     d[j] <- dj 
    }
	if( ncol(y) == 2 ){
		 n <- rowSums(y) 
		 i <- y1 == 0 
         j <- y1 == n 
         Y <- y1/n
         d <- y1*log(y1/(n*yhat)) + (n-y1)*log((1-Y)/(1-yhat)) 
	     di <- -n[i]*log(1-yhat[i])
	     dj <- -n[j]*log(yhat[j])
	     d[i] <- di
	     d[j] <- dj
	}
    
    vu <- yhat*(1-yhat)/n      
    pear.res <- (Y-yhat) / sqrt( vu )     
    Xp <- sum(pear.res^2)
    
    
    w <- diag( n*yhat / (1-yhat) )
    inv <- MASS::ginv(t(X) %*% w %*% X)
    H <- w %*% X %*% inv %*% t(X) 
    
    lev <- diag(H)
    cook <-  lev * ( pear.res/(1-lev) )^2 / ncol(X)
    
   fit <- list(coefficients = beta, 
                 se = SE , 
                 vcov = out$vcov, 
                 vcov0=vcov0, 
                 #vcov1=vcov1,
                 fitted.values = yhat, 
                 df = nrow(X) - ncol(X), 
                 loglik= -out$value,
                 deviance = 2*sum( d ),  # = 2*out$value when n_i = 1, all i 
                 dev.resid = sign(Y-yhat)*sqrt(2*d),
                 residuals = pear.res,
                 X2 = Xp,
                 hat.matrix = H,
                 cook.distance = cook,   
                 convergence = out$convergence,
                 barrier.value = out$barrier.value,
                 outer.iterations = out$outer.iterations, 
                 Active = active, 
                 formula = formula,
                 tol = tol, 
                 start.beta = start.beta)
    
    class(fit) <- "lbreg"
    return(fit)
 
}

