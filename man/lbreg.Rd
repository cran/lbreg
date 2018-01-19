\name{lbreg}
\alias{lbreg}
\title{Log-Binomial regression}

\description{Fitting a Log-Binomial Regression Model}

\usage{
lbreg(formula, data, start.beta, tol=0.9999, delta=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{formula}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
\item{data}{ 
an optional data frame containing the variables in the model. 
If not found in data, the variables are taken from environment(formula), typically the environment from which lbreg is called.
}
\item{start.beta}{
starting values for the parameters in the linear predictor. If missing, the default value explained in 
Andrade and Andrade (2018) is used according to the choice of \code{delta}.
}
\item{tol}{
defaults to 0.9999; threshold for declaring a probability on the boundary (p = 1). 
}
\item{delta}{
defaults to 1. See reference below.
}
\item{\dots}{ not used. }

}
 

\details{
This function uses \code{\link{constrOptim}} with the BFGS method in order
to perform maximum likelihood estimation of the log-binomial regression model
as described in the reference below. When the MLE is the interior of the parameter space results 
should agree with \code{glm(...,family=binomial(link='log'))}.
\code{lbreg} uses the adaptive logarithimic barrier algorithm rather than 
iteratively weighted least squares (\code{glm}).
}
\value{
\item{Active}{matrix of active constraints.}
\item{barrier.value}{same as in \code{constrOptim}.}
\item{coefficients}{named vector of estimated regression coefficients.}
\item{convergence}{same as in \code{constrOptim}.}
\item{call}{the matched call.}
\item{cook.distance}{Cook's distance.} 
\item{data}{the data argument.}
\item{deviance}{residual deviance.}
\item{dev.resid}{deviance residuals.}
\item{fitted.values}{fitted probabilities.}
\item{formula}{the formula supplied.}
\item{hat.matrix}{hat matrix for GLMs (whose diagonal contains leverage values).}
\item{loglik}{maximized loglikelihood.}
\item{outer.iterations}{same as in \code{constrOptim}.}
\item{residuals}{Pearson residuals.}
\item{se}{standard errors of estimated coefficients.}
\item{start.beta}{starting values used by \code{constrOptim}.}
\item{vcov}{variance-covariance matrix of estimates.}
\item{vcov0}{inverse of observed Fisher information; should be equal to vcov if 
there are no active constraints (Active = NULL).}
\item{X2}{sum of squared residuals (variance-inflation estimate (dispersion) = X2/df).}
}
\references{
Andrade, BB; Andrade JML (forthcoming) Some results for Maximum Likelihood Estimation of Adjusted Relative Risks. 
Communications in Statistics - Theory and Methods.
}
\author{
\packageAuthor{lbreg}
}
%\note{
%%  ~~further notes~~
%}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{glm}} \code{(family=binomial(link='log'))}, \code{\link{relrisk}}
}

\examples{
require(lbreg)

# data preparation
data(PCS)  # ungrouped data
w <- PCS
w <- w[,-1]
w$race <- factor(w$race)
w$dpros <- factor(w$dpros)
w$dcaps <- factor(w$dcaps)

# log-binomial regression
fm <- lbreg(tumor ~ ., data=w)
fm
coef(fm)
summary(fm)


# grouped data
require(lbreg)
data(Caesarian)
m1 <- lbreg( cbind(n1, n0) ~ RISK + NPLAN + ANTIB, data=Caesarian)
summary(m1)

# dispersion estimate based on deviance residuals
sum(m1$dev.res^2)  
# dispersion estimate based on Pearson residuals (reported in the summary above)
sum(m1$residuals^2)/(8-4)  

predict(m1, newdata=data.frame(RISK=0, NPLAN=1, ANTIB=1))

# m0 <- glm( cbind(n1, n0) ~ RISK + NPLAN + ANTIB, data=Dat, family=binomial(link='log'))
# summary(m0)


}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
