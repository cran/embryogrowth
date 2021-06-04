#' Result of the fit using the nest database with linear progression
#' @title Fit using the nest database
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @docType data
#' @name resultNest_4p_SSM_Linear
#' @encoding UTF-8
#' @description Fit using the nest database
#' @keywords datasets
#' @usage resultNest_4p_SSM_Linear
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' pfixed <- NULL
#' M0 = 0
#' ############################################################################
#' # 4 parameters SSM
#' ############################################################################
#' x <- c('DHA' = 64.868697530424186, 'DHH' = 673.18292743646771, 
#'        'T12H' = 400.90952554047749, 'Rho25' = 82.217237723502123)
#' resultNest_4p_SSM_Linear <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, integral=integral.linear, M0=M0, 
#' 	hatchling.metric=c(Mean=39.33, SD=1.92)/39.33)
#' }
#' @format A list with fitted information about data(nest)
NULL
