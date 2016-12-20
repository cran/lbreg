summary.lbreg <- function(object, ...)
{
  se <- object$se
  zval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               z.value = zval,
               p.value = 2*pnorm(-abs(zval),mean=0,sd=1))
  res <- list(call=object$call, coefficients=TAB
              #AIC=object$AIC,
              #Deviance=object$deviance
  )
  class(res) <- "summary.lbreg"
  res
}
