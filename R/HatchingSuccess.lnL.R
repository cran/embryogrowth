#' HatchingSuccess.lnL return -log likelihood of the data and the parameters
#' @title Return -log likelihood of the data and the parameters
#' @author Marc Girondot
#' @return Return -log likelihood of the data and the parameters
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
#'                                !is.na(Total) & Total != 0 & 
#'                                !is.na(NotHatched) & !is.na(Hatched))
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
#' t <- seq(from=20, to=40, by=0.1)
#' CIq <- predict(g, temperature=t)
#' 
#' par(mar=c(4, 4, 1, 1), +0.4)
#' plot(g)
#' }
#' @export



HatchingSuccess.lnL <- function(par, data, fixed.parameters=NULL, 
                  column.Incubation.temperature="Incubation.temperature",
                  column.Hatched="Hatched",
                  column.NotHatched="NotHatched") {
  # S.low, P.low logistique basse
  # S.high, deltaP logistic haute
  # MaxHS
  model <- HatchingSuccess.model(par=c(par, fixed.parameters), temperature=data[, column.Incubation.temperature])
  model <- ifelse(model == 0, 1E-9, model)
  model <- ifelse(model == 1, 1 - 1E-9, model)
  lnL <- -sum(dbinom(x = data[, column.Hatched], size=data[, column.Hatched]+data[, column.NotHatched], prob=model, log=TRUE))
  if (is.infinite(lnL)) {
    # print(par)
    lnL <- 1E9
  }
  return(lnL)
}

