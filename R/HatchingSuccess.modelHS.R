#' modelHS returns the hatching success according the set of parameters and temperatures
#' @title Return the hatching success according the set of parameters and temperatures
#' @author Marc Girondot
#' @return Return the hatching success according the set of parameters and temperatures
#' @param x A set of parameters.
#' @param temperature A vector of temperatures.
#' @description Set of functions to study the hatching success.\cr
#' @family Hatching success
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' totalIncubation_Cc <- subset(DatabaseTSD, 
#'                              Species=="Caretta caretta" & 
#'                                Note != "Sinusoidal pattern" & 
#'                                !is.na(Total) & Total != 0)
#' 
#' par <- c(S.low=0.5, S.high=0.3, 
#'          P.low=25, deltaP=10, MaxHS=logit(0.8))
#'          
#' lnLHS(x=par, data=totalIncubation_Cc)
#' 
#' g <- fitHS(par=par, data=totalIncubation_Cc)
#' 
#' lnLHS(x=g$par, data=totalIncubation_Cc)
#' 
#' t <- seq(from=20, to=40, by=0.1)
#' CIq <- predict(g, temperature=t)
#' 
#' par(mar=c(4, 4, 1, 1), +0.4)
#' plot(x=t, 
#'      y=modelHS(x=g$par, temperature = t), 
#'      bty="n", las=1, type="n", ylim=c(0,1), xlim=c(20, 40), 
#'      xlab="Incubation temperature", ylab="Hatching success")
#' with(totalIncubation_Cc, points(x = Incubation.temperature, 
#'      y = Hatched/(Hatched+NotHatched), 
#'      col="red", pch=19))
#' lines(x = t, y=CIq["2.5%", ], lty=2)
#' lines(x = t, y=CIq["97.5%", ], lty=2)
#' lines(x = t, y=CIq["50%", ], lty=1)
#' }
#' @export



modelHS <- function(x, temperature) {
  P.high <- x["P.low"]+abs(x["deltaP"])
  model <- (1/(1+exp((1/abs(x["S.low"]))*(x["P.low"]-temperature))))*(1/(1+exp((1/-abs(x["S.high"]))*(P.high-temperature))))*invlogit(x["MaxHS"])
  return(model)
}


