#' Result of the fit using the nest database
#' @title Fit using the nest database
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name resultNest_6p_SSM
#' @encoding UTF-8
#' @description Fit using the nest database
#' @references Girondot, M., Monsinjon, J., Guillon, J.-M., Submitted. 
#'             Delimitation of the embryonic thermosensitive period for 
#'             sex determination using an embryo growth model reveals a 
#'             potential bias for sex ratio prediction in turtles.
#' @keywords datasets
#' @usage resultNest_6p_SSM
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' x <- structure(c(104.954347370542, 3447.10062406071, 661.269363920423, 
#'  96.3871849546537, 306.456389026151, 232.105840347154), .Names = c("DHA", 
#'  "DHH", "DHL", "DT", "T12L", "Rho25"))
#' pfixed <- c(rK=1.208968)
#' resultNest_6p_SSM <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=0.3470893, 
#' 	test=c(Mean=39.33, SD=1.92))
#' plotR(result=resultNest_6p_SSM, show.hist = TRUE,
#'              ylim=c(0, 8), curves="ML")
#' }
#' @format A list with fitted information about data(nest)
NULL