#' Timeseries of constant temperatures for nests
#' @title Timeseries of constant temperatures for nests
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name tempConst
#' @description Timeseries of temperatures for nests
#' @references Girondot, M., Kaska, Y., 2014. A model to predict the thermal reaction 
#' norm for the embryo growth rate from field data. Journal of Thermal Biology 45, 96-102.
#' @keywords datasets
#' @usage tempConst
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' # Same as:
#' # GenerateConstInc(durations = rep(104*60*24, 11),
#' # temperatures = 25:35,
#' # names = paste0("T",25:35))
#' data(tempConst)
#' tempConst_f <- FormatNests(tempConst)
#' 
#' data(nest)
#' formated <- FormatNests(nest)
#' x <- structure(c(109.683413821537, 614.969219372661, 306.386903812694, 
#'  229.003478775323), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#'  
#'  # See the stages dataset examples for justification of M0 and rK
#'  
#' pfixed <- c(rK=1.208968)
#' resultNest_4p_SSM <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=0.3470893, 
#' 	test=c(Mean=39.33, SD=1.92))
#' 	
#' plotR(result=resultNest_4p_SSM, show.hist = TRUE,
#'              ylim=c(0, 8), curves="ML quantiles")
#' 
#' # Now use the fited parameters from resultNest_4p_SSM with  
#' # the constant incubation temperatures:
#' 
#' plot(resultNest_4p_SSM, temperatures=tempConst_f,  
#' 	stopattest=TRUE, series="T30", xlim=c(0,50),  
#' 	ylimT=c(22, 32), test=c(Mean=39.33, SD=1.92), 
#' 	embryo.stages="Caretta caretta.SCL")
#' 	
#' plot(resultNest_4p_SSM, temperatures=tempConst_f,  
#' 	stopattest=TRUE, series="T25", xlim=c(0,120),  
#' 	ylimT=c(22, 32), test=c(Mean=39.33, SD=1.92), 
#' 	embryo.stages="Caretta caretta.SCL")
#' 	
#' }
#' @format A dataframe with raw data.
NULL
