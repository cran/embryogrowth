#' Result of the mcmc using the nest database for epsilon parameter
#' @title Result of the mcmc using the nest database
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name resultNest_mcmc_4p_epsilon
#' @encoding UTF-8
#' @description Fit using the nest database
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage resultNest_mcmc_4p_epsilon
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
#' resultNest_4p_SSM4p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p_SSM4p)
#' resultNest_4p_epsilon <- resultNest_4p_SSM4p
#' resultNest_4p_epsilon$fixed.parameters <- c(resultNest_4p_epsilon$par, 
#'                    resultNest_4p_epsilon$fixed.parameters)
#' resultNest_4p_epsilon$par <- c(epsilon = 0)
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p_epsilon, accept = TRUE)
#' resultNest_mcmc_4p_epsilon <- GRTRN_MHmcmc(result = resultNest_4p_epsilon, 
#'            n.iter = 10000, parametersMCMC = pMCMC,
#'            n.chains = 1, n.adapt = 0, thin = 1, trace = TRUE, parallel = TRUE)
#' data(resultNest_mcmc_4p_epsilon)
#' plot(resultNest_mcmc_4p_epsilon, parameters="epsilon", xlim=c(-11, 11), las=1)
#' plotR(resultNest_4p_epsilon, SE=c(epsilon = unname(resultNest_mcmc_4p_epsilon$SD)), 
#'            ylim=c(0, 3), las=1)
#' }
#' @format A list of class mcmcComposite with mcmc result for data(nest) with 4 parameters and Gompertz model of growth
NULL
