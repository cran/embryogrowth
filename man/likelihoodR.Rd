% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihoodR.R
\name{likelihoodR}
\alias{likelihoodR}
\title{Estimate the likelihood of a set of parameters for nest incubation data}
\usage{
likelihoodR(
  result = NULL,
  parameters = NULL,
  fixed.parameters = NULL,
  temperatures = NULL,
  integral = NULL,
  derivate = NULL,
  hatchling.metric = NULL,
  M0 = NULL,
  hessian = FALSE,
  weight = NULL,
  parallel = TRUE,
  echo = TRUE
)
}
\arguments{
\item{result}{A object obtained after searchR or likelihoodR}

\item{parameters}{A set of parameters}

\item{fixed.parameters}{A set of parameters that will not be changed}

\item{temperatures}{Timeseries of temperatures}

\item{integral}{Function used to fit embryo growth: integral.Gompertz, integral.exponential or integral.linear}

\item{derivate}{Function used to fit embryo growth: dydt.Gompertz, dydt.exponential or dydt.linear. It will replace the one in NestsResult.}

\item{hatchling.metric}{Mean and SD of size of hatchlings}

\item{M0}{Measure of hatchling size or mass proxi at laying date}

\item{hessian}{If TRUE, the hessian matrix is estimated and the SE of parameters estimated.}

\item{weight}{A named vector of the weight for each nest for likelihood estimation}

\item{parallel}{If true, try to use several cores using parallel computing.}

\item{echo}{If FALSE, does not display the result.}
}
\value{
A result object
}
\description{
Estimate the likelihood of a set of parameters for nest incubation data
}
\details{
likelihoodR estimates the likelihood of a set of parameters for nest incubation data
}
\examples{
\dontrun{
library(embryogrowth)
data(nest)
formated <- FormatNests(nest)
# The initial parameters value can be:
# "T12H", "DHA",  "DHH", "Rho25"
# Or
# "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
# K for Gompertz must be set as fixed parameter or being a constant K  
# or relative to the hatchling size rK
x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
# pfixed <- c(K=82.33) or rK=82.33/39.33
pfixed <- c(rK=2.093313)
# K or rK are not used for integral.linear or integral.exponential
LresultNest_4p <- likelihoodR(parameters=x, fixed.parameters=pfixed,  
	temperatures=formated, integral=integral.Gompertz, M0=1.7,  
	hatchling.metric=c(Mean=39.33, SD=1.92))
data(resultNest_4p_SSM)
LresultNest_4p <- likelihoodR(result=resultNest_4p_SSM)
}
}
\author{
Marc Girondot \email{marc.girondot@gmail.com}
}
