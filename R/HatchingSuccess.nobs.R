#' nobs.NestsResult Return number of observations of a fit
#' @title Return number of observations of a fit
#' @author Marc Girondot
#' @return Return number of observations of a fit
#' @param object The return of a fit done with fitHS.
#' @param ... Not used
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
#' HatchingSuccesss.lnL(par=par, data=totalIncubation_Cc)
#' 
#' g <- HatchingSuccesss.fit(par=par, data=totalIncubation_Cc)
#' 
#' HatchingSuccesss.lnL(par=g$par, data=totalIncubation_Cc)
#' 
#' plot(g)
#' }
#' @method nobs HatchingSuccess
#' @export

nobs.HatchingSuccess <- function(object, ...) {
  return(nrow(object$data))
}

