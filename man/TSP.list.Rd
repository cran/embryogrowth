% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TSP.list.R
\docType{data}
\name{TSP.list}
\alias{TSP.list}
\title{Database of thermosensitive period of development for sex determination}
\format{
A list with dataframes including attributes
}
\usage{
TSP.list
}
\description{
Database of thermosensitive period of development for sex determination.\cr
This database can be used with the functions plot() or info.nests().\cr
The attributes TSP.begin.stages and TSP.end.stages for each dataframe give
respectively the first and the last stages for TSP. Then the metrics for the limits of TSP
are the average sizes before and after the TSP (see example, below).\cr
If the metric for the stages before the TSP or after the TSP is not known, it will use the
available information.
}
\details{
Database of thermosensitive period of development for sex determination
}
\examples{
\dontrun{
library(embryogrowth)
data(TSP.list)
names(TSP.list)
reference <- "Emys_orbicularis.mass"
metric <- TSP.list[[reference]]
TSP.begin <- attributes(TSP.list[[reference]])$TSP.begin.stages
TSP.end <- attributes(TSP.list[[reference]])$TSP.end.stages
# Metric at the beginning of the TSP
del <- ifelse(all(metric$stages == TSP.begin - 1)==FALSE, 0, 1)
(metric$metric[metric$stages == TSP.begin - del] +
    metric$metric[metric$stages == TSP.begin]) / 2
# Metric at the end of the TSP
del <- ifelse(all(metric$stages == TSP.begin + 1)==FALSE, 0, 1)
(metric$metric[metric$stages == TSP.end] +
        metric$metric[metric$stages == del + TSP.end]) / 2
}
}
\references{
\insertRef{1056}{embryogrowth}\cr
\insertRef{10870}{embryogrowth}\cr
\insertRef{10620}{embryogrowth}\cr
}
\seealso{
Other Functions for temperature-dependent sex determination: 
\code{\link{DatabaseTSD}},
\code{\link{DatabaseTSD.version}()},
\code{\link{P_TRT}()},
\code{\link{ROSIE}},
\code{\link{ROSIE.version}()},
\code{\link{plot.tsd}()},
\code{\link{predict.tsd}()},
\code{\link{stages}},
\code{\link{tsd}()},
\code{\link{tsd_MHmcmc}()},
\code{\link{tsd_MHmcmc_p}()}
}
\author{
Marc Girondot \email{marc.girondot@universite-paris-saclay.fr}
}
\concept{Functions for temperature-dependent sex determination}
\keyword{datasets}
