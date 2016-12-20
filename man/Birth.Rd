\name{Birth}
\alias{Birth}
\docType{data}
\title{Birth Weight Data}
\description{Data used by Wacholder (1986) to illustrate the use of log binomial 
regression for estimating adjusted relative risks  of  a  low-birthweight  baby. }
\usage{data("Birth")}
\format{
  A data frame with 900 observations on the following 5 variables.
  \describe{
\item{\code{lowbw}}{low birth weight delivery (1=yes) }
\item{\code{alc}}{mother's alcohol drinking frequency (1=Light, 2=Moderate, 3=Heavy)}
\item{\code{smo}}{mother smoked (1=no)}		
\item{\code{soc}}{mother's social status (1=I and II (lower), 2=III (middle), 3=IV and V (upper))}
  }
}
%%\details{}
\source{
Stata's online manual  http://www.stata.com/manuals13/rbinreg.pdf
}
\references{
Wright JT, Waterson EJ, Barrison PJ, et al. (1983). Alcohol consumption, pregnancy and	low birthweight. Lancet 1:663-665.
}
\examples{
data(Birth)
dim(Birth)
names(Birth)

}
\keyword{datasets}
