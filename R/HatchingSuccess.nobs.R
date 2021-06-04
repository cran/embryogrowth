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
#'          P.low=25, deltaP=10, MaxHS=0.8)
#'          
#' HatchingSuccess.lnL(par=par, data=totalIncubation_Cc)
#' 
#' g <- HatchingSuccess.fit(par=par, data=totalIncubation_Cc)
#' 
#' HatchingSuccess.lnL(par=g$par, data=totalIncubation_Cc)
#' 
#' plot(g)
#' 
#' nobs(g)
#' }
#' @method nobs HatchingSuccess
#' @export

nobs.HatchingSuccess <- function(object, ...) {
  return(sum(object$data[, object$column.Hatched] + object$data[, object$column.NotHatched]))
}

