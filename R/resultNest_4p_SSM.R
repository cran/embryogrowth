#' Result of the fit using the nest database
#' @title Fit using the nest database
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @docType data
#' @name resultNest_4p_SSM
#' @encoding UTF-8
#' @description Fit using the nest database
#' @references Girondot M, Monsinjon J, Guillon J-M (2018) Delimitation of the 
#' embryonic thermosensitive period for sex determination using an embryo 
#' growth model reveals a potential bias for sex ratio prediction in turtles. 
#' Journal of Thermal Biology 73: 32-40
#' @keywords datasets
#' @usage resultNest_4p_SSM
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' x <- structure(c(109.683413821537, 614.969219372661, 306.386903812694, 
#'  229.003478775323), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' pfixed <- c(rK=1.208968)
#' resultNest_4p_SSM <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, integral=integral.Gompertz, M0=0.3470893, 
#' 	hatchling.metric=c(Mean=39.33, SD=1.92))
#' plotR(result=resultNest_4p_SSM, show.hist = TRUE,
#'              ylim=c(0, 8), curve="ML quantiles")
#' }
#' @format A list with fitted information about data(nest)
NULL
