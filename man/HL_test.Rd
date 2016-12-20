\name{HL_test}
\alias{HL_test}

\title{Hosmer-Lemeshow Goodness of Fit Test}
\description{
The HL decile-of-risk test.  Validity of the test
assumes that the number of covariate patterns
is close to the number of observations which is violated 
when many observations have the same covariate pattern. 
and several ties will impact the required ordering and grouping
of observations. 
This is less likely when there is at least one 
continuous covariate. 
}
\usage{ HL_test(object, g = 10) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{  object of class 'lbreg'. }
  \item{g}{  number of groups }
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A list with elements
\item{X2}{HL statistic}
\item{pvalue}{p-value for the test from Chi Squared with df = g-2}
}
\references{
Hosmer D W, Lemeshow S 2000. Applied Logistic Regression. New York, USA: John Wiley and Sons.
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

# data preparation
data(PCS)
w <- PCS
w <- w[,-1]
w$race <- factor(w$race)
w$dpros <- factor(w$dpros)
w$dcaps <- factor(w$dcaps)

fm <- lbreg(tumor ~ ., data=w)

HL_test(fm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
