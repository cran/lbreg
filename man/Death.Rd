\name{Death}
\alias{Death}
\docType{data}
\title{Death Penalty Data}
\description{
See references. 
}
\usage{data("Death")}
\format{
  A data frame with 147 observations on the following 6 variables.
  \describe{
    \item{\code{death}}{death = 1, life in prison = 0}
    \item{\code{blackd}}{black defendant = 1}
    \item{\code{whitvic}}{white victim = 1}
    \item{\code{serious}}{a measure of crime seriousness}
    \item{\code{culp}}{a measure of culpability}
    \item{\code{serious2}}{another measure of crime seriousness}
  }
}

%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
SAS Institute Inc. (2006). Logistic regression using the SAS system: Theory and application.
SAS Publishing, Cary, NC: SAS Institute Inc; http://ftp.sas.com/~samples/A55770
}
\references{
Petersen MR, Deddens JA (2010). Maximum Likelihood Estimation of the Log-Binomial
Model. Communications in Statistics: Theory and Methods, 39, 874-883.
}
\examples{
data(Death)
dim(Death)
names(Death)
}
\keyword{datasets}
