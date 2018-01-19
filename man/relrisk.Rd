\name{relrisk}
\alias{relrisk}

\title{Regression Adjusted Relative Risks
%%  ~~function to do ... ~~
}
\description{
This function calculates the relative risks RR adjusted for 
covariates (acting on a previous log-binomial regression fit) 
and confidence intervals (by default 95 percent) for the estimated RR. 
The confidence interval is calculated from the log(RR) and backtransformed. 
}
\usage{relrisk(object, alpha = 0.05, dispersion = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class 'lbreg'.}
  \item{alpha}{1 - desired confidence level.}
  \item{dispersion}{logical. \code{TRUE} if standard errors should be adjusted for dispersion estimate based on Pearson residuals.} 
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{value}{table with estimated relative risks, lower and upper bounds of condifidence intervals.}
}
\references{
Andrade, BB; Andrade JML (forthcoming) Some results for Maximum Likelihood Estimation of Adjusted Relative Risks. 
Communications in Statistics - Theory and Methods.
}
\author{
\packageAuthor{lbreg}}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lbreg}}
}
\examples{
require(lbreg)

# ungrouped data
# data preparation
data(PCS)
w <- PCS
w <- w[,-1]
w$race <- factor(w$race)
w$dpros <- factor(w$dpros)
w$dcaps <- factor(w$dcaps)

# log-binomial regression
fm <- lbreg(tumor ~ ., data=w)

# relative risks
relrisk(fm)
relrisk(fm, alpha=.10)

# grouped data 
require(lbreg)
data(Caesarian)
m1 <- lbreg( cbind(n1, n0) ~ RISK + NPLAN + ANTIB, data=Caesarian)
relrisk(m1)
relrisk(m1, dispersion=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
