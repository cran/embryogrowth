#' Result of the fit using the nest database with anchored parameters
#' @title Fit using the nest database with anchored parameters
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name resultNest_newp
#' @description Fit using the nest database with anchored parameters
#' @references Girondot, M. & Kaska, Y. In press. A model to predict the thermal 
#'          reaction norm for the embryo growth rate from field data. Journal of
#'          Thermal Biology. In press
#' @keywords datasets
#' @usage resultNest_newp
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' newp <- GenerateAnchor(nests=formated, number.anchors=7)
#' pfixed <- c(rK=2.093313)
#' resultNest_newp <- searchR(parameters=newp, fixed.parameters=pfixed, 
#'   temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92), method = "BFGS", maxiter = 200)
#' data(resultNest_newp)
#' }
#' @format A list with fitted information about data(nest) with anchored parameters
NULL
