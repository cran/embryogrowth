#' predict.HatchingSuccess returns prediction based on a model fitted with fitHS
#' @title Return prediction based on a model fitted with fitHS
#' @author Marc Girondot
#' @return Return a matrix with prediction based on a model fitted with fitHS
#' @param temperature A vector of temperatures.
#' @param probs Quantiles.
#' @param replicates Number of replicates to estimate the confidence interval.
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
#' @method predict HatchingSuccess
#' @export




predict.HatchingSuccess <- function(object, ..., 
                                    temperature=NULL, 
                                    probs=c(0.025, 0.5, 0.975), 
                                    replicates=1000) {
  SE <- object$SE
  par <- object$par
  
  # mathessian <- object$hessian
  if (is.null(replicates)) replicates <- 1
  if (is.null(temperature)) temperature <- object$data$Incubation.temperature

  CI <- matrix(data = NA , ncol=replicates, nrow=length(temperature))
  
  CI[, 1] <- modelHS(par, temperature)
  if (replicates >1) {
  for (c in 2:replicates) {
    x <- structure(rnorm(n = 5, mean=par, sd=SE), .Names=names(par))
    CI[, c] <- modelHS(x, temperature)
  }
  }
  
  CIq <- apply(CI, MARGIN = 1, FUN = function(x) {quantile(x, probs = probs)})
  if (class(CIq)=="numeric") CIq <- matrix(data=CIq, ncol=length(temperature))
  colnames(CIq) <- as.character(temperature)
  
  return(CIq)
}

