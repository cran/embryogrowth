% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nobs.NestsResult.R
\name{nobs.NestsResult}
\alias{nobs.NestsResult}
\title{Return number of observations of a fit}
\usage{
\method{nobs}{NestsResult}(object, ...)
}
\arguments{
\item{object}{A result file generated by searchR}

\item{...}{Not used}
}
\value{
Return number of observations of a fit
}
\description{
Return number of observations of a fit.\cr
This function is used for bbmle::ICtb().
}
\details{
nobs.NestsResult Return number of observations of a fit
}
\examples{
\dontrun{
library(embryogrowth)
data(resultNest_4p_SSM)
logLik(resultNest_4p_SSM)
AIC(resultNest_4p_SSM)
nobs(resultNest_4p_SSM)
}
}
\author{
Marc Girondot
}
