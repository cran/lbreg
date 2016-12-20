relrisk <-
function(object, alpha=0.05)
{
  rr <- exp(object$coefficients)
  zcrit <- qnorm(1 - alpha/2)
  rr.ci.e.low <- exp( object$coefficients - zcrit * object$se )   # CI of the form exp[ b -/+ z SE(b) ]
  rr.ci.e.up <- exp( object$coefficients + zcrit * object$se )   # CI of the form exp[ b -/+ z SE(b) ]
  TAB <- cbind( RR = rr,
                lowCI = rr.ci.e.low,
                upCI = rr.ci.e.up  )
  cat("CI limits correspond to level", 1-alpha, "\n") 
  return(TAB)
}
