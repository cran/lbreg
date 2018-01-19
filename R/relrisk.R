relrisk <- function(object, alpha=0.05, dispersion=FALSE)
{
  if(dispersion){
	   m <- sqrt(object$X2/object$df)	
  	   cat("(Standard Errors have been multiplied by sqrt of dispersion estimate.)")
	   cat("\n")
   }else{
	   m <- 1
   }
  
  rr <- exp(object$coefficients)
  zcrit <- qnorm(1 - alpha/2)
  
  # CI of the form exp[ b -/+ z SE(b) ]
  rr.ci.e.low <- exp( object$coefficients - zcrit * m * object$se )   
  rr.ci.e.up <- exp( object$coefficients + zcrit * m * object$se ) 
  TAB <- cbind( RR = rr,
                lowCI = rr.ci.e.low,
                upCI = rr.ci.e.up  )
  cat("CI limits correspond to level", 1-alpha, "\n") 
  return(TAB)
  }
