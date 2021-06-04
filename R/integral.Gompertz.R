#' integral.Gompertz returns the derivative of the Gompertz function.
#' @title Return the result of the Gompertz function
#' @author Marc Girondot
#' @return A list with the derivative
#' @param t The time in any unit
#' @param size The current size
#' @param parms A vector with alpha and K values being c(alpha=x1, K=x2)
#' @description Return the result of the Gompertz function as a data.frame with two columns, time and metric\cr
#' integral.Gompertz(t, size, parms)
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
#' # K or rK are not used for dydt.linear or dydt.exponential
#' resultNest_4p_SSM <- searchR(parameters=x, fixed.parameters=pfixed,  
#' 	temperatures=formated, integral=integral.Gompertz, M0=1.7,  
#' 	hatchling.metric=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p_SSM)
#' }
#' @export

integral.Gompertz <- function(t, size, parms) {
		return(data.frame(time=t, metric=parms["K"]*exp(log(size/parms["K"]) * exp(-parms["alpha"] * t))))
		}
