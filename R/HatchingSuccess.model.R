#' HatchingSuccess.model returns the hatching success according the set of parameters and temperatures
#' @title Return the hatching success according the set of parameters and temperatures
#' @author Marc Girondot
#' @return Return the hatching success according the set of parameters and temperatures
#' @param par A set of parameters.
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
#' HatchingSuccess.lnL(par=par, data=totalIncubation_Cc)
#' 
#' g <- HatchingSuccess.fit(par=par, data=totalIncubation_Cc)
#' 
#' HatchingSuccess.lnL(par=g$par, data=totalIncubation_Cc)
#' 
#' plot(g)
#' }
#' @export



HatchingSuccess.model <- function(par, temperature) {
  P.high <- par["P.low"]+abs(par["deltaP"])
  model <- (1/(1+exp((1/abs(par["S.low"]))*(par["P.low"]-temperature))))*(1/(1+exp((1/-abs(par["S.high"]))*(P.high-temperature))))*invlogit(par["MaxHS"])
  return(model)
}


