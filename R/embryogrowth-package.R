#' Fit a parametric function that describes dependency of embryo growth to temperature
#'
#' \tabular{ll}{
#'  Package: \tab embryogrowth\cr
#'  Type: \tab Package\cr
#'  Version: \tab 5.0 - build 396\cr
#'  Date: \tab 2014-08-15\cr
#'  License: \tab GPL (>= 2)\cr
#'  LazyLoad: \tab yes\cr
#'  }
#' @title The package embryogrowth
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType package
#' @name embryogrowth-package
#' @description A package to fit growth of embryos
#' @references Girondot, M. & Kaska, Y. In press. A model to predict the thermal 
#'          reaction norm for the embryo growth rate from field data. Journal of
#'          Thermal Biology. In press
#' @seealso Delmas, V., Prevot-Julliard, A.-C., Pieau, C. & Girondot, M. 2008. 
#'          A mechanistic model of temperature-dependent sex determination 
#'          in a Chelonian, the European pond turtle. Functional 
#'          Ecology, 22, 84-93.
#' @seealso Girondot, M., Ben Hassine, S., Sellos, C., Godfrey, M. & Guillon, 
#'          J.-M. 2010. Modeling thermal influence on animal growth and sex 
#'          determination in Reptiles: being closer of the target gives new 
#'          views. Sexual Development, 4, 29-38.
#' @seealso Girondot, M. 1999. Statistical description of temperature-dependent 
#'          sex determination using maximum likelihood. Evolutionary Ecology 
#'          Research, 1, 479-486.
#' @keywords Temperature Embryo Ecology Growth Gompertz
#' @examples
#' \dontrun{
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(115.758929130522, 428.649022170996, 503.687251738993, 
#' 12.2621455821612, 306.308841227278, 116.35048615105), .Names = c("DHA", 
#' "DHH", "DHL", "DT", "T12L", "Rho25"))
#' # or
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
#' out <- as.mcmc(result_mcmc_4p)
#' # This out obtained after as.mcmc can be used with coda package
#' # plot() can use the direct output of embryogrowth_MHmcmc() function.
#' plot(result_mcmc_4p, parameters=1, xlim=c(0,550))
#' plot(result_mcmc_4p, parameters=3, xlim=c(290,320))
#' # summary() permits to get rapidly the standard errors for parameters
#' summary(result_mcmc_4p)
#' # The batch standard error procedure is usually thought to  
#' # be not as accurate as the time series methods.
#' se <- result_mcmc_4p$BatchSE
#' # or
#' se <- result_mcmc_4p$TimeSeriesSE
#' }

NULL
