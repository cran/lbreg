summary.lbreg <- function(object, ...)
{
  se <- object$se
  zval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               z.value = zval,
               p.value = 2*pnorm(-abs(zval),mean=0,sd=1))
  res <- list(call=object$call, 
              coefficients=TAB, 
              loglik=object$loglik, 
              df=object$df,
              deviance=object$deviance,   
              residuals=object$residuals,
              X2 = object$X2,
              yclass = attr(object$terms, 'dataClasses')[1]  # identifies grouped data
              #AIC=object$AIC,
  )
  class(res) <- "summary.lbreg"
  res
}
