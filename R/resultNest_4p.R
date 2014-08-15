#' Result of the fit using the nest database
#' @title Fit using the nest database
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name resultNest_4p
#' @description Fit using the nest database
#' @references Girondot, M. & Kaska, Y. Submitted. A model to predict 
#'             temperature dependency on embryo growth rate and incubation
#'             duration from field data.
#' @keywords datasets
#' @usage resultNest_4p
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' resultNest_4p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92), method = "BFGS", maxiter = 200)
#' }
#' @format A list with fitted information about data(nest)
NULL
