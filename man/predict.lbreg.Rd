\name{predict.lbreg}
\alias{predict.lbreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict method for Log-Binomial regression.
}
\description{
Predicted values based on 'lbreg' object.
}
\usage{
\method{predict}{lbreg}(object, newdata, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class inheriting from "lbreg"
}
  \item{newdata}{  
  a data frame with covariate values with which to predict. If omitted, the fitted probabilities are returned.
}
  \item{\dots}{ not used
%%     ~~Describe \code{\dots} here~~
}
}
\details{
If newdata is omitted the predictions are simply the fitted values stored in the object supplied.
}
\value{
\item{Active}{active restrictions (taking newdata into account).}
\item{coef.pred}{regression coefficients re-estimated to satisfy possibly new restrictions imposed by newdata. See reference below.}
\item{convergence}{same as in the object supplied.}
\item{se.pred}{estimated standard errors of predictions.}
\item{tol}{same as in the object supplied.}
\item{ypred}{predicted probabilities for newdata.}
}
\references{
Andrade, BB; Andrade JML (2018) Some results for Maximum Likelihood Estimation of Adjusted Relative Risks. 
Communications in Statistics - Theory and Methods.

}
\author{
\packageAuthor{lbreg}
}
%\note{
%%  ~~further notes~~
%}

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
novo <- data.frame(age=c(41, 32), race=c(1,2), dpros=c(2,4), 
                   dcaps=c(1,1), psa=c(7.24,3.25), vol=c(4.3,5.6), 
                   gleason=c(2,8))
predict(fm, newdata=novo)
}
