#' logLik.HatchingSuccess returns -log L of a fit
#' @title Return -log L of a fit
#' @author Marc Girondot
#' @return Return -log L of a fit
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
#' g <- HatchingSuccess.fit(par=par, data=totalIncubation_Cc)
#' 
#' HatchingSuccesss.lnL(par=g$par, data=totalIncubation_Cc)
#' 
#' t <- seq(from=20, to=40, by=0.1)
#' CIq <- predict(g, temperature=t)
#' 
#' par(mar=c(4, 4, 1, 1), +0.4)
#' plot(g)
#' }
#' @method logLik HatchingSuccess
#' @export



logLik.HatchingSuccess <- function(object, ...) {
  L <- -object$value
  attributes(L) <- list(df=length(object$par), 
                        nobs=nrow(object$data))
  return(L)
}

