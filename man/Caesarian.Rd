\name{Caesarian}
\alias{Caesarian}
\docType{data}
\title{Caesarian Infection Dataset}
\description{Adapted dataset from Fahrmeir et al (2013): grouped data on infections of 251 mothers after a 
C-section collected at the clinical center of the University of Munich.}
\usage{data("Caesarian")}
\format{
  A data frame with 7 rows and 5 variables.
  \describe{
    \item{\code{n1}}{Caesarians with infections.}
    \item{\code{n0}}{Caesarians without infections.}
    \item{\code{NPLAN}}{ = 1 if C-section was not planned.}
    \item{\code{RISK}}{= 1 if risk factors existed.}
    \item{\code{ANTIB}}{ = 1 if antibiotics were administered as prophylaxis.}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
http://www.uni-goettingen.de/de/551625.html
}
\references{
Fahrmeir, L., Kneib, Th., Lang, S., Marx, B. (2013) Regression - Models, Methods and Applications. Springer. 
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Caesarian)
Caesarian
# no observations for case (RISK=0, NPLAN=1, ANTIB=1)
y = Caesarian[,1:2]
cbind(Caesarian[,3:5], total=rowSums(y))
colSums(y)
}
\keyword{datasets}
