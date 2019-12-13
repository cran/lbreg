\name{lbreg-package}
\alias{lbreg-package}
\docType{package}
\title{
\packageTitle{lbreg}
}
\description{
\packageDescription{lbreg}

Package \code{lbreg} performs maximum likelihood estimation of Log-Binomial Regression. 
The main functions are \code{\link{lbreg}} which provides a shortcut to \code{\link{constrOptim}} 
to estimate LBR coefficients and \code{\link{relrisk}} which takes lbreg results to 
produce estimated relative risks and associated confidence intervals and prediction. Results
differ from \code{\link{glm}} when the MLE is on the boundary of the parameter space as 
explained in the reference below (Andrade, Andrade (2018)).


}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{lbreg}
\packageIndices{lbreg}

}
\author{
\packageAuthor{lbreg}

Maintainer: \packageMaintainer{lbreg}
}
\references{
Andrade, BB; Andrade JML (2018) Some results for Maximum Likelihood Estimation of Adjusted Relative Risks. 
Communications in Statistics - Theory and Methods.

}



