\name{Evans}
\alias{Evans}
\docType{data}
\title{Evans County dataset}
\description{Data from cohort study in which white males in Evans County 
were followed for 7 years, with coronary
heart disease as the outcome of interest.}
\usage{data("Evans")}
\format{
  A data frame with 609 observations on the following 9 variables.
  \describe{
\item{\code{CDH}}{outcome variable; 1 = coronary heart disease}
\item{\code{CAT}}{1 = high, 0 = normal catecholamine level}
\item{\code{AGE}}{age (in years)}
\item{\code{CHL}}{cholesterol, mg/dl}
\item{\code{SMK}}{1 = subject has ever smoked}
\item{\code{ECG}}{1 = presence of electrocardiogram abnormality}
\item{\code{DBP}}{diastolic blood pressure, mmHg}
\item{\code{SBP}}{systolic blood pressure, mmHg}
\item{\code{HPT}}{1 = SBP greater than or equal to 160 or DBP greater than or equal to 95}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
http://web1.sph.emory.edu/dkleinb/logreg3.htm#data
}
\references{
D. Kleinbaum and M. Klein (2010) Survival Analysis: A Self-Learning Text. 3rd ed. Springer. 
}
\examples{
data(Evans)
dim(Evans)
names(Evans)

}
\keyword{datasets}
