\name{Heart}
\alias{Heart}
\docType{data}
\title{Heart Dataset}
\description{Heart attack data from the ASSENT-2 study.}
\usage{data("Heart")}
\format{
  A data frame with 16,949 observations on the following 5 variables.
  \describe{
    \item{\code{Heart}}{binary response; 1 = death}
    \item{\code{age}}{categorized into <65, 65-75 or >75 years}
    \item{\code{severity}}{Killip class I, II, or III/IV}
    \item{\code{region}}{code for three USA regions}
    \item{\code{onset}}{treatment delay categorized into <2, 2-4 or >4 hours}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
http://biostatistics.oxfordjournals.org/content/13/1/179/suppl/DC1
}
\references{
ASSENT-2 INVESTIGATORS (1999). Single-bolus tenecteplase compared with front-loaded alteplase in acute
myocardial infarction: the ASSENT-2 double-blind randomised trial. Lancet 354, 716-722.


Ian C. Marschner and Alexandra C. Gillett (2012)
Relative risk regression: reliable and flexible methods for log-binomial models.
Biostatistics 13, 179-192 
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Heart)
dim(Heart)
names(Heart)

}
\keyword{datasets}
