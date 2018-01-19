print.summary.lbreg <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")

  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)

  cat("\n")
  cat("Residual deviance: ", x$deviance," on", x$df, " degrees of freedom")
  cat("\n")
  
  if(x$yclass == "nmatrix.2"){
  cat("Dispersion Estimate: ", x$X2/x$df)
  cat("\n")
  cat("(Std. Errors have *not* been multiplied by sqrt of dispersion estimate.)")
  cat("\n")
  }
}
