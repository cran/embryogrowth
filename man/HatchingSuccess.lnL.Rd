% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HatchingSuccess.lnL.R
\name{HatchingSuccess.lnL}
\alias{HatchingSuccess.lnL}
\title{Return -log likelihood of the data and the parameters}
\usage{
HatchingSuccess.lnL(
  par,
  data,
  fixed.parameters = NULL,
  column.Incubation.temperature = "Incubation.temperature",
  column.Hatched = "Hatched",
  column.NotHatched = "NotHatched"
)
}
\arguments{
\item{par}{A set of parameters.}

\item{data}{A dataset in a data.frame with a least three columns: Incubation.temperature, Hatched and NotHatched}

\item{fixed.parameters}{A set of parameters that must not be fitted.}

\item{column.Incubation.temperature}{Name of the column with incubation temperatures}

\item{column.Hatched}{Name of the column with hatched number}

\item{column.NotHatched}{Name of the column with not hatched number}
}
\value{
Return -log likelihood of the data and the parameters
}
\description{
Set of functions to study the hatching success.\cr
}
\details{
HatchingSuccess.lnL return -log likelihood of the data and the parameters
}
\examples{
\dontrun{
library(embryogrowth)
totalIncubation_Cc <- subset(DatabaseTSD, 
                             Species=="Caretta caretta" & 
                               Note != "Sinusoidal pattern" & 
                               !is.na(Total) & Total != 0 & 
                               !is.na(NotHatched) & !is.na(Hatched))

par <- c(S.low=0.5, S.high=0.3, 
         P.low=25, deltaP=10, MaxHS=0.8)
         
HatchingSuccess.lnL(par=par, data=totalIncubation_Cc)

g <- HatchingSuccess.fit(par=par, data=totalIncubation_Cc)

HatchingSuccess.lnL(par=g$par, data=totalIncubation_Cc)

t <- seq(from=20, to=40, by=0.1)
CIq <- predict(g, temperature=t)

par(mar=c(4, 4, 1, 1), +0.4)
plot(g)
}
}
\seealso{
Other Hatching success: 
\code{\link{HatchingSuccess.MHmcmc}()},
\code{\link{HatchingSuccess.MHmcmc_p}()},
\code{\link{HatchingSuccess.fit}()},
\code{\link{HatchingSuccess.model}()},
\code{\link{logLik.HatchingSuccess}()},
\code{\link{nobs.HatchingSuccess}()},
\code{\link{predict.HatchingSuccess}()}
}
\author{
Marc Girondot
}
\concept{Hatching success}
