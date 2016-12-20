\name{PCS}
\alias{PCS}
\docType{data}
\title{PCS Dataset}
\description{Prostate Cancer Study}
\usage{data("PCS")}
\format{
  A data frame with 380 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{Identification Code; 1 - 380}
    \item{\code{tumor}}{Tumor Penetration of Prostatic Capsule,	0 = No Penetration}
    \item{\code{age}}{in years}
    \item{\code{race}}{Race; 1= White, 2 = Black}
    \item{\code{dpros}}{Results of the Digital Rectal Exam, 4 levels}
    \item{\code{dcaps}}{Detection of Capsular Involvement in Rectal Exam; 	1 = No, 2 = Yes	}
    \item{\code{psa}}{antigen mg/ml}
    \item{\code{vol}}{Tumor Volume Obtained from Ultrasound, cm3}
    \item{\code{gleason}}{Total Gleason Score; 0 - 10}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
https://www.umass.edu/statdata/statdata/data/pros.txt
}
\references{
Hosmer and Lemeshow (2000) Applied Logistic Regression, Wiley.  

%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(PCS)
## View(PCS)
## str(PCS) ; plot(PCS) ...
}
\keyword{datasets}
