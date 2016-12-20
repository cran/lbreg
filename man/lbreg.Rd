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
Andrade and Andrade (submitted) is used according to the choice of \code{delta}.
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
This function is basically a wrapper to \code{\link{constrOptim}} with BFGS method in order
to perform maximum likelihood estimation of the lob-binomial regression model
as described in the reference below. When the MLE is the interior of the parameter space results should agree with 
\code{glm(formula, family=binomial(link='log'), data)}.
\code{lbreg} uses the adaptive logarithimic barrier algorithm rather than 
iteratively weighted least squares (\code{glm}).
This package is work in progress and it currently provides less functionality than \code{\link{glm}}.
}
\value{
\item{Active}{matrix of active constraints.}
\item{barrier.value}{same as in \code{constrOptim}.}
\item{coefficients}{named vector of estimated regression coefficients.}
\item{convergence}{same as in \code{constrOptim}.}
\item{call}{the matched call.}
\item{data}{the data argument.}
\item{fitted.values}{fitted probabilities.}
\item{formula}{the formula supplied.}
\item{loglik}{maximized loglikelihood.}
\item{outer.iterations}{same as in \code{constrOptim}.}
\item{residuals}{raw residuals (yi - pi).}
\item{se}{standard errors of estimates.}
\item{start.beta}{starting values supplied.}
\item{vcov}{var-cov matrix of estimates.}
\item{vcov0}{inverse of observed Fisher information; should be equal to vcov if 
there are no active constraints (Active = NULL).}
}
\references{
Andrade, BB; Andrade JML (submitted) Some results for Maximum Likelihood Estimation of Adjusted Relative Risks. 
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
data(PCS)
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


}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
