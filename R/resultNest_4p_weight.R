#' Result of the fit using the nest database with weight
#' @title Fit using the nest database with weight
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name resultNest_4p_weight
#' @description Fit using the nest database with weight
#' @references Girondot, M. & Kaska, Y. In press. A model to predict the thermal 
#'          reaction norm for the embryo growth rate from field data. Journal of
#'          Thermal Biology. In press
#' @keywords datasets
#' @usage resultNest_4p_weight
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' w <- weightmaxentropy(formated, control_plot=list(xlim=c(20,36)))
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' # K or rK are not used for dydt.linear or dydt.exponential
#' resultNest_4p_weight <- searchR(parameters=x,  
#' 	fixed.parameters=pfixed, temperatures=formated,  
#' 	derivate=dydt.Gompertz, M0=1.7, test=c(Mean=39.33, SD=1.92),  
#' 	method = "BFGS", maxiter = 200, weight=w)
#' data(resultNest_4p_weight)
#' plotR(resultNest_4p_weight, ylim=c(0,0.50), xlim=c(15, 35))
#' }
#' @format A list with fitted information about data(nest)
NULL
