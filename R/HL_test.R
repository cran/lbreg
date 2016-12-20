
HL_test <- function( object, g=10){
	
   if(class(object) != "lbreg") warning("HL.test was meant for _lbreg_ objects; also OK for some _glm's_")
   
   DATA <- object$data
   y <- DATA[,as.character(formula(object)[[2]])]
   yhat <- fitted(object) 

   br <- unique( quantile(yhat, probs=seq(from=0, to=1,  by=1/g) ) )
   qy <- cut(yhat, breaks=br, labels=FALSE)

   Obs0 <- xtabs(1 - y ~ qy)
   Exp0 <- xtabs(1 - yhat ~ qy)
  
   Obs1 <- xtabs(y ~ qy)
   Exp1 <- xtabs(yhat ~ qy)
  
   hl <- sum( (Obs0-Exp0)^2 / Exp0  + (Obs1-Exp1)^2 / Exp1 )

   return( list(X2=hl, pvalue= 1-pchisq(hl, df=g-2), g=g) )
}
