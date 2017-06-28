# Function to perform LBR via constrained optimization
# (basically a wrapper to 'constrOptim' which is already part of R base) 

lbregEst <- function(x, y, start.beta, tol, ...)
{
    X <- x
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
    out$fitted.values <- yhat
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
    Wtil = diag(p) - proj(t(A))
    W = Wtil[1:p,1:(p-r)]
    if( max( A%*%W ) > sqrt( 1e-6 ) ) { warning("max AW is > 1e-6") }
    meat = try( MASS::ginv(t(W) %*% J %*% W) )
      if(class(meat) != "try-error"){
         new.vcov <- W %*% meat %*% t(W)
         out$vcov <- new.vcov
      }else{
         warning("could not invert (W'JW) --- 'vcov' ignores active constraints")
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
     
      
    return( list(coefficients = beta, 
                 se = SE , 
                 vcov = out$vcov, 
                 vcov0=vcov0, 
                 #vcov1=vcov1,
                 fitted.values = yhat, 
                 df = nrow(X) - ncol(X), 
                 loglik= -out$value,
                 convergence = out$convergence,
                 barrier.value = out$barrier.value,
                 outer.iterations = out$outer.iterations, 
                 Active = active, 
                 formula = formula,
                 tol = tol, 
                 start.beta = start.beta)
     ) 
}

