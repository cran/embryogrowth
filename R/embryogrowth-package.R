#' Fit a parametric function that describes dependency of embryo growth to temperature
#'
#' \tabular{ll}{
#'  Package: \tab embryogrowth\cr
#'  Type: \tab Package\cr
#'  Version: \tab 10.2 build 2094\cr
#'  Date: \tab 2025-06-16\cr
#'  License: \tab GPL (>= 2)\cr
#'  LazyLoad: \tab yes\cr
#'  }
#' @title The package embryogrowth
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @name embryogrowth-package
#' @description Tools to analyze the embryo growth and the sexualisation thermal reaction norms.\cr
#' The latest version of this package can always been installed using:\cr
#' install.packages("http://marc.girondot.free.fr/CRAN/HelpersMG.tar.gz", repos=NULL, type="source")\cr
#' install.packages("http://marc.girondot.free.fr/CRAN/embryogrowth.tar.gz", repos=NULL, type="source")\cr
#' \if{html}{\figure{E.png}{options: alt="embryogrowth logo"}}
#' \if{latex}{\figure{E.png}}
#' @references
#' \insertRef{1515}{embryogrowth}\cr
#' \insertRef{3534}{embryogrowth}\cr
#' \insertRef{9039}{embryogrowth}\cr
#' \insertRef{8589}{embryogrowth}\cr
#' \insertRef{10871}{embryogrowth}\cr
#' \insertRef{8566}{embryogrowth}\cr
#' \insertRef{11124}{embryogrowth}\cr
#' \insertRef{10620}{embryogrowth}\cr
#' \insertRef{11754}{embryogrowth}\cr
#' \insertRef{12168}{embryogrowth}\cr
#' \insertRef{10870}{embryogrowth}\cr
#' \insertRef{13669}{embryogrowth}\cr
#' \insertRef{5790}{embryogrowth}\cr
#' \insertRef{13271}{embryogrowth}\cr
#' \insertRef{11893}{embryogrowth}\cr
#' @importFrom Rdpack reprompt
#' @keywords Temperature Embryo Ecology Growth Gompertz Sex-determination
#' @examples
#' \dontrun{
#' library("embryogrowth")
#' packageVersion("embryogrowth")
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
#' x <- structure(c(118.431040984352, 498.205702157603, 306.056280989839, 
#' 118.189669472381), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' 
#' ################################################################################
#' #
#' # The values of rK=2.093313 and M0=1.7 were used in 
#' # Girondot, M. & Kaska, Y. 2014. A model to predict the thermal 
#' # reaction norm for the embryo growth rate from field data. Journal of
#' # Thermal Biology. 45, 96-102.
#' #
#' # Based on recent analysis on table of development for both Emys orbicularis and 
#' # Caretta caretta, best value for rK should be 1.209 and M0 should be 0.34.
#' # Girondot M, Monsinjon J, Guillon J-M (2018) Delimitation of the embryonic 
#' # thermosensitive period for sex determination using an embryo growth model 
#' # reveals a potential bias for sex ratio prediction in turtles. Journal of 
#' # Thermal Biology 73: 32-40 
#' #
#' # See the example in the stages datasets
#' # 
#' ################################################################################
#' 
#' resultNest_4p_SSM <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, integral=integral.Gompertz, M0=1.7, 
#' 	hatchling.metric=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p_SSM)
#' par(mar=c(4, 4, 1, 1))
#' plot(resultNest_4p_SSM$data, bty="n", las=1, 
#'      xlab="Days of incubation", ylab="Temperatures in °C", 
#'      series="all", 
#'      type="l", xlim=c(0,70),ylim=c(20, 35))
#' par(mar=c(4, 4, 1, 1))
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p_SSM, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' resultNest_mcmc_4p_SSM <- GRTRN_MHmcmc(result=resultNest_4p_SSM,  
#' 	parametersMCMC=pMCMC, n.iter=10000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(resultNest_mcmc_4p_SSM)
#' out <- as.mcmc(resultNest_mcmc_4p_SSM)
#' # This out obtained after as.mcmc can be used with coda package
#' # plot() can use the direct output of GRTRN_MHmcmc() function.
#' plot(resultNest_mcmc_4p_SSM, parameters=1, xlim=c(0,550))
#' plot(resultNest_mcmc_4p_SSM, parameters=3, xlim=c(290,320))
#' # But rather than to use the SD for each parameter independantly, it is
#' # more logical to estimate the distribution of the curves
#' new_result <- ChangeSSM(resultmcmc = resultNest_mcmc_4p_SSM, result = resultNest_4p_SSM,
#'                         temperatures = seq(from = 20, to = 35, by = 0.1), 
#'                         initial.parameters = NULL)
#' par(mar=c(4, 4, 1, 5)+0.4)
#' 
#' plotR(result = resultNest_4p_SSM, parameters = new_result$par, 
#'            ylabH = "Temperatures\ndensity", ylimH=c(0, 0.3), atH=c(0, 0.1, 0.2), 
#'            ylim=c(0, 3), show.hist=TRUE)
#'       
#' # Beautiful density plots
#' 
#' plotR(result = resultNest_4p_SSM, 
#'              resultmcmc=resultNest_mcmc_4p_SSM, 
#'              ylim=c(0, 8), 
#'              curve = "MCMC quantiles", show.density=TRUE)
#' 
#' plotR(resultNest_6p_SSM, resultmcmc=resultNest_mcmc_6p_SSM, 
#'             ylim=c(0, 8), show.density=TRUE, show.hist=TRUE, 
#'             curve = "MCMC quantiles", 
#'             ylimH=c(0,0.5), atH=c(0, 0.1, 0.2))
#' 
#' # How many times this package has been download
#' library(cranlogs)
#' embryogrowth <- cran_downloads("embryogrowth", from = "2014-08-16", 
#'                             to = Sys.Date() - 1) 
#' sum(embryogrowth$count)
#' plot(embryogrowth$date, embryogrowth$count, type="l", bty="n")
#' }

NULL
