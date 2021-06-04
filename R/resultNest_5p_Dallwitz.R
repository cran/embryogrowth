#' Result of the fit using the nest database
#' @title Fit using the nest database
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @docType data
#' @name resultNest_5p_Dallwitz
#' @encoding UTF-8
#' @description Fit using the nest database
#' @references Girondot M, Monsinjon J, Guillon J-M (2018) Delimitation of the embryonic thermosensitive period for sex determination using an embryo growth model reveals a potential bias for sex ratio prediction in turtles. Journal of Thermal Biology 73: 32-40 
#' @keywords datasets
#' @usage resultNest_5p_Dallwitz
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' x <- structure(c(4.91191231405918, 12.7453211281394, 31.2670410811077, 
#'  5.7449376569153, -0.825689964543813), .Names = c("Dallwitz_b1",
#'  "Dallwitz_b2", "Dallwitz_b3", "Dallwitz_b4", "Dallwitz_b5"))
#' pfixed <- c(rK=1.208968)
#' resultNest_5p_Dallwitz <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, integral=integral.Gompertz, M0=0.3470893, 
#' 	hatchling.metric=c(Mean=39.33, SD=1.92))
#' plotR(result=resultNest_5p_Dallwitz, show.hist = TRUE,
#'              ylim=c(0, 8), curves="ML quantiles")
#' }
#' @format A list with fitted information about data(nest)
NULL
