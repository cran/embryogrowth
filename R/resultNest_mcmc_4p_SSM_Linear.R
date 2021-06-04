#' Result of the mcmc using the nest database
#' @title Result of the mcmc using the nest database
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @docType data
#' @name resultNest_mcmc_4p_SSM_Linear
#' @encoding UTF-8
#' @description Fit using the nest database
#' @keywords datasets
#' @usage resultNest_mcmc_4p_SSM_Linear
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
#' data(resultNest_4p_SSM_Linear)
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p_SSM_Linear, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' resultNest_mcmc_4p_SSM_Linear <- GRTRN_MHmcmc(result=resultNest_4p_SSM_Linear, 
#'  adaptive = TRUE, 
#' 	parametersMCMC=pMCMC, n.iter=10000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(resultNest_mcmc_4p_SSM_Linear)
#' 1-rejectionRate(as.mcmc(resultNest_mcmc_4p_SSM_Linear))
#' as.parameters(resultNest_mcmc_4p_SSM_Linear)
#' layout(mat=matrix(1:4, nrow = 2))
#' plot(resultNest_mcmc_4p_SSM_Linear, parameters = "all", scale.prior = TRUE, las = 1)
#' layout(mat=1)
#' plotR(resultNest_4p_SSM_Linear, resultmcmc=resultNest_mcmc_4p_SSM_Linear, ylim=c(0,4), 
#'          main="Schoolfield, Sharpe & Magnuson 4-parameters", show.density=TRUE)
#' }
#' @format A list of class mcmcComposite with mcmc result for data(nest) with 4 parameters and linear model of growth
NULL
