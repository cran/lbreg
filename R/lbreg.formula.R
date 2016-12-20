lbreg.formula <- function(formula, data=list(), ...)
{
  data <- na.omit(data)
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  
  est <- lbreg.default(x, y, ...)

  est$call <- match.call()
  est$formula <- formula
  est$data <- data
  est
}
