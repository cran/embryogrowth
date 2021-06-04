#' Result of the mcmc using the nest database
#' @title Result of the mcmc using the nest database
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @docType data
#' @name resultNest_mcmc_4p_SSM
#' @encoding UTF-8
#' @description Fit using the nest database
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage resultNest_mcmc_4p_SSM
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(118.431040984352, 498.205702157603, 306.056280989839, 
#' 118.189669472381), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' resultNest_4p_SSM <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, integral=integral.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p_SSM)
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p_SSM, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' resultNest_mcmc_4p_SSM <- GRTRN_MHmcmc(result=resultNest_4p_SSM, 
#'  adaptive = TRUE, 
#' 	parametersMCMC=pMCMC, n.iter=10000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(resultNest_mcmc_4p_SSM)
#' 1-rejectionRate(as.mcmc(resultNest_mcmc_4p_SSM))
#' as.parameters(resultNest_mcmc_4p_SSM)
#' layout(mat=matrix(1:4, nrow = 2))
#' plot(resultNest_mcmc_4p_SSM, parameters = "all", scale.prior = TRUE, las = 1)
#' layout(mat=1)
#' plotR(resultNest_4p_SSM, resultmcmc=resultNest_mcmc_4p_SSM, ylim=c(0,4), 
#'          main="Schoolfield, Sharpe & Magnuson 4-parameters", show.density=TRUE)
#' }
#' @format A list of class mcmcComposite with mcmc result for data(nest) with 4 parameters and Gompertz model of growth
NULL
