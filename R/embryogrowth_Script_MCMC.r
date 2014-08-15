#' embryogrowth_MHmcmc runs the Metropolis-Hastings algorithm for data (Bayesian MCMC)
#' @title Run the Metropolis-Hastings algorithm for data
#' @author Marc Girondot
#' @return A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
#' @param n.iter Number of iterations for each step
#' @param parametersMCMC A set of parameters used as initial point for searching with information on priors
#' @param result An object obtained after a SearchR fit
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param trace True or False, shows progress
#' @param batchSize Number of observations to include in each batch fo SE estimation
#' @param parallel If true, try to use several cores using parallel computing.
#' @description Run the Metropolis-Hastings algorithm for data.\cr
#' Deeply modified from a MCMC script by Olivier Martin (INRA, Paris-Grignon).\cr
#' The number of iterations is n.iter+n.adapt+1 because the initial likelihood is also displayed.\cr
#' I recommend that thin=1 because the method to estimate SE uses resampling.\cr
#' If initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
#' Note that resultLnL has all the likelihoods, not only those defined by n.adapt and thin.\cr
#' To get the SE from result_mcmc <- embryogrowth_MHmcmc(result=try), use:\cr
#' result_mcmc$BatchSE or result_mcmc$TimeSeriesSE\cr
#' The batch standard error procedure is usually thought to be not as accurate as the time series methods.\cr
#' Based on Jones, Haran, Caffo and Neath (2005), the batch size should be equal to sqrt(n.iter).\cr
#' Jones, G.L., Haran, M., Caffo, B.S. and Neath, R. (2006) Fixed Width Output Analysis for Markov chain Monte Carlo , Journal of the American Statistical Association, 101:1537-1547.\cr
#' coda package is necessary for this function.\cr
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "T12H", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' resultNest_4p <- searchR(parameters=x, fixed.parameters=pfixed,  
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7,  
#' 	test=c(Mean=39.33, SD=1.92), method = "BFGS", maxiter = 200)
#' data(resultNest_4p)
#' pMCMC <- embryogrowth_MHmcmc_p(resultNest_4p, accept=TRUE)
#' # Take care, it can be very long; several days
#' result_mcmc_4p <- embryogrowth_MHmcmc(result=resultNest_4p, 
#' 		parametersMCMC=pMCMC, n.iter=10000, n.chains = 1,  
#' 		n.adapt = 0, thin=1, trace=TRUE)
#' data(result_mcmc_4p)
#' out <- as.mcmc(result_mcmc_4p)
#' # This out can be used with coda package
#' # plot() can use the direct output of embryogrowth_MHmcmc() function.
#' plot(result_mcmc_4p, parameters=1, xlim=c(0,550))
#' plot(result_mcmc_4p, parameters=3, xlim=c(290,320))
#' # summary() permits to get rapidly the standard errors for parameters
#' summary(result_mcmc_4p)
#' # They are store in the result also. Two SE are estimated using or 
#' # batch method or time-series SE:
#' # The batch standard error procedure is usually thought to be not 
#' # as accurate as the time series methods.
#' se1 <- result_mcmc_4p$BatchSE
#' se2 <- result_mcmc_4p$TimeSeriesSE
#' }
#' @export

embryogrowth_MHmcmc <- function(result=stop("An output from searchR must be provided"), n.iter=10000, 
parametersMCMC=stop("Priors must be given. Use embryogrowth_MHmcmc_p()"), n.chains = 1, n.adapt = 0, 
thin=1, trace=FALSE, batchSize=sqrt(n.iter), parallel=TRUE)
{

  if (any(installed.packages()[,1]=="coda")) {
    require("coda") } else {
      warning("coda package is necessary for this function")
      return()
    }

parameters <- parametersMCMC

print(parameters)

# 29/1/2014; Ajout de result$weight
out <- .MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, 
n.chains = n.chains, n.adapt = n.adapt, thin=thin, trace=trace, 
	data=list(temperatures=result$data, derivate=result$derivate, 
		testuse=result$test, M0=result$M0, fixed.parameters=result$fixed.parameters, 
		parallel=parallel, weight=result$weight), likelihood=get(".fonctionMCMC"))

if (batchSize>=n.iter/2) {
  print("batchSize cannot be larger than half the number of iterations.")
  rese <- rep(NA, dim(parametersMCMC)[1])
  names(rese) <- rownames(parametersMCMC)
  out <- c(out, SE=list(rese))
} else {
  out <- c(out, BatchSE=list(batchSE(out$resultMCMC, batchSize=batchSize)))
}

class(out) <- "mcmcComposite"

out <- c(out, TimeSeriesSE=list(summary(out)$statistics[,4]))

class(out) <- "mcmcComposite"

return(out)
}
