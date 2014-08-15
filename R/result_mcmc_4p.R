#' Result of the mcmc using the nest database
#' @title Result of the mcmc using the nest database
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name result_mcmc_4p
#' @description Fit using the nest database
#' @references Girondot, M. & Kaska, Y. Submitted. A model to predict 
#'             temperature dependency on embryo growth rate and incubation
#'             duration from field data.
#' @keywords datasets
#' @usage result_mcmc_4p
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
#' data(resultNest_4p)
#' pMCMC <- embryogrowth_MHmcmc_p(resultNest_4p, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' result_mcmc_4p <- embryogrowth_MHmcmc(result=resultNest_4p,  
#' 	parametersMCMC=pMCMC, n.iter=10000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(result_mcmc_4p)
#' plot(result_mcmc_4p, parameters="T12H", main="", xlim=c(290, 320), bty="n")
#' plotR(resultNest_4p, SE=result_mcmc_4p$TimeSeriesSE, ylim=c(0,0.3), las=1)
#' }
#' @format A list of class mcmcComposite with mcmc result for data(nest) with 4 parameters and Gompertz model of growth
NULL
