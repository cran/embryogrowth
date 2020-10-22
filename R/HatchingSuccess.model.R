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
#'          P.low=25, deltaP=10, MaxHS=0.8)
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
  P.low <- par["P.low"]
  deltaP <- abs(par["deltaP"])
  P.high <- abs(par["P.high"])
  
  if (is.na(P.low)) P.low <- P.high - deltaP
  if (is.na(P.high)) P.high <- P.low + deltaP
  if (is.na(deltaP)) deltaP <- P.high - P.low
  
  if (is.na(par["K1.low"])) par <- c(par, K1.low=0)
  if (is.na(par["K2.low"])) par <- c(par, K2.low=0)
  if (is.na(par["K1.high"])) par <- c(par, K1.high=0)
  if (is.na(par["K2.high"])) par <- c(par, K2.high=0)
  
  m.low <- flexit(x=temperature, 
               P=P.low, S=par["S.low"], 
               K1=par["K1.low"], K2=par["K2.low"])
  m.high <- flexit(x=temperature, 
               P=P.high, S=-par["S.high"], 
               K1=par["K1.high"], K2=par["K2.high"])
  
  
  model <- m.low*m.high*par["MaxHS"]
  return(model)
}


