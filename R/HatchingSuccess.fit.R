#' HatchingSuccess.fit fits a hatching success model to data
#' @title Fit a hatching success model to data using maximum likelihood
#' @author Marc Girondot
#' @return Return a object of class HatchingSuccess
#' @param par A set of parameters.
#' @param data A dataset in a data.frame with a least three columns: Incubation.temperature, Hatched and NotHatched
#' @param fixed.parameters A set of parameters that must not be fitted.
#' @param column.Incubation.temperature Name of the column with incubation temperatures
#' @param column.Hatched Name of the column with hatched number
#' @param column.NotHatched Name of the column with not hatched number
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
#' par(mar=c(4, 4, 1, 1), +0.4)
#' plot(g)
#' }
#' @export


HatchingSuccess.fit <- function(par, data, fixed.parameters=NULL, 
                  column.Incubation.temperature="Incubation.temperature",
                  column.Hatched="Hatched",
                  column.NotHatched="NotHatched") {
  

  g <- suppressWarnings(optim(par=par, fn=HatchingSuccess.lnL, 
             fixed.parameters=fixed.parameters, 
             column.Incubation.temperature=column.Incubation.temperature,
             column.Hatched=column.Hatched,
             column.NotHatched=column.NotHatched,
             data=data, hessian = FALSE))
  
  g$par["S.low"] <- abs(g$par["S.low"])
  g$par["S.high"] <- abs(g$par["S.high"])

  g <- suppressWarnings(optim(par=g$par, fn=HatchingSuccess.lnL, 
                              fixed.parameters=fixed.parameters, 
                              column.Incubation.temperature=column.Incubation.temperature,
                              column.Hatched=column.Hatched,
                              column.NotHatched=column.NotHatched,
                              data=data, hessian = TRUE))
  
  g$fixed.parameters <- fixed.parameters
  g$SE <- SEfromHessian(g$hessian)
  g$AIC <- 2*g$value+2*length(g$par)
  g$data <- data
  g$column.Incubation.temperature <- column.Incubation.temperature
  g$column.Hatched <- column.Hatched
  g$column.NotHatched <- column.NotHatched
    
  class(g) <- "HatchingSuccess"
  return(g)
}


