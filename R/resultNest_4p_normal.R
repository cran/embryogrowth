#' Result of the fit using the nest database using asymmetric normal function
#' @title Result of the fit using the nest database using asymmetric normal function
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @docType data
#' @name resultNest_4p_normal
#' @encoding UTF-8
#' @description Fit using the nest database using asymmetric normal function
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage resultNest_4p_normal
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' x <- ChangeSSM(temperatures = (200:350)/10,
#'                parameters = resultNest_4p_SSM$par,
#'                initial.parameters = structure(c(3, 7, 11, 32), 
#'                                .Names = c("Scale", "sdL", "sdH", "Peak")), 
#'                control=list(maxit=1000))
#' pfixed <- c(rK=2.093313)
#' resultNest_4p_normal <- searchR(parameters=x$par, fixed.parameters=pfixed, 
#'                          temperatures=formated, integral=integral.Gompertz, M0=1.7, 
#'                          hatchling.metric=c(Mean=39.33, SD=1.92))
#' plotR(resultNest_4p_normal, ylim=c(0, 3))
#' plotR(resultNest_4p_normal, ylim=c(0, 3), ylimH = c(0, 0.9), show.hist=TRUE)
#' compare_AIC(SSM=resultNest_4p_SSM, Asymmetric.normal=resultNest_4p_normal)
#' }
#' @format A list with fitted information about data(nest)
NULL
