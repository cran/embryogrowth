#' Result of the mcmc using the nest database with anchored parameters
#' @title Result of the mcmc using the nest database with anchored parameters
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name result_mcmc_newp
#' @description Fit using the nest database with anchored parameters
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage result_mcmc_newp
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # the package polynom must be present
#' if (is.element('polynom', installed.packages()[,1]) == FALSE) {
#' install.packages('polynom')
#' }
#' newp <- GenerateAnchor(nests=formated, number.anchors=7)
#' pfixed <- c(rK=2.093313)
#' resultNest_newp <- searchR(parameters=newp, fixed.parameters=pfixed, 
#'   temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#'   test=c(Mean=39.33, SD=1.92), method = "BFGS", maxiter = 200)
#' data(resultNest_newp)
#' pMCMC <- TRN_MHmcmc_p(resultNest_newp, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' result_mcmc_newp <- GRTRN_MHmcmc(result=resultNest_newp,  
#'   parametersMCMC=pMCMC, n.iter=1000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(result_mcmc_newp)
#' data(resultNest_4p)
#' newp <- GenerateAnchor(nests=resultNest_4p, number.anchors=7)
#' # Here the confidence interval is built based on anchored parameters
#' plotR_hist(resultNest_4p, parameters=newp, SE=result_mcmc_newp$TimeSeriesSE, 
#'  ylim=c(0,0.4), ylimH=c(0,0.4))
#' # Here the confidence interval is built based on parametric SSM equation
#' data(result_mcmc_4p)
#' plotR_hist(resultNest_4p, SE=result_mcmc_4p$TimeSeriesSE,
#'  ylim=c(0,0.4), ylimH=c(0,0.4))
#' plot(result_mcmc_newp, las=1, xlim=c(0,30), parameters="294")
#' plot(result_mcmc_newp, las=1, xlim=c(0,30), parameters="296.333333333333")
#' plot(result_mcmc_newp, las=1, xlim=c(0,30), parameters=3)
#' }
#' @format A list of class mcmcComposite with mcmc result for data(nest) with anchored parameters
NULL
