lbreg.formula <- function(formula, data=list(), ...)
{
  data <- na.omit(data)
  mf <- model.frame(formula=formula, data=data)
  terms <- attr(mf, "terms")
    
  fit <- lbreg.default(x=model.matrix(terms, mf), y=model.response(mf), ...)  
  
  fit$terms <- terms
  fit$call <- match.call()
  fit$formula <- formula
  fit$data <- data
  fit
}
