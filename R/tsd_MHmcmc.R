#' tsd_MHmcmc runs the Metropolis-Hastings algorithm for tsd (Bayesian MCMC)
#' @title Run the Metropolis-Hastings algorithm for tsd
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
#' @description Run the Metropolis-Hastings algorithm for tsd.\cr
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
#' eo <- subset(STSRE_TSD, Species=="Emys orbicularis", c("Males", "Females", 
#'                                        "Incubation.temperature"))
#' eo_logistic <- tsd(eo)
#' pMCMC <- tsd_MHmcmc_p(eo_logistic, accept=TRUE)
#' # Take care, it can be very long; several days
#' result_mcmc_tsd <- tsd_MHmcmc(result=eo_logistic, 
#' 		parametersMCMC=pMCMC, n.iter=10000, n.chains = 1,  
#' 		n.adapt = 0, thin=1, trace=TRUE)
#' # summary() permits to get rapidly the standard errors for parameters
#' summary(result_mcmc_tsd)
#' # They are store in the result also. Two SE are estimated using or 
#' # batch method or time-series SE:
#' # The batch standard error procedure is usually thought to be not 
#' # as accurate as the time series methods.
#' se1 <- result_mcmc_tsd$BatchSE
#' se2 <- result_mcmc_tsd$TimeSeriesSE
#' plot(result_mcmc_tsd, parameters="S", scale.prior=TRUE, xlim=c(-3, 3), las=1)
#' plot(result_mcmc_tsd, parameters="P", scale.prior=TRUE, xlim=c(25, 35), las=1)
#' plot(eo_logistic, se=se2)
#' }
#' @export

tsd_MHmcmc <- function(result=stop("An output from tsd() must be provided"), n.iter=10000, 
parametersMCMC=stop("Priors must be given. Use tsd_MHmcmc_p()"), n.chains = 1, n.adapt = 0, 
thin=1, trace=FALSE, batchSize=sqrt(n.iter)) {

# result=eo_logistic; parametersMCMC=pMCMC; n.iter=10000; n.chains = 1;  n.adapt = 0; thin=1; trace=TRUE; batchSize=sqrt(n.iter)
  
  require("coda")
  
parameters <- parametersMCMC

print(parameters)

# 29/1/2014; Ajout de result$weight
out <- .MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, 
n.chains = n.chains, n.adapt = n.adapt, thin=thin, trace=trace, 
	data=list(males=result$males, N=result$N, temperatures=result$temperatures, 
            equation=result$equation), likelihood=get(".fonctiontsdMCMC"))

if (batchSize>=n.iter/2) {
  print("batchSize cannot be larger than half the number of iterations.")
  rese <- rep(NA, dim(parametersMCMC)[1])
  names(rese) <- rownames(parametersMCMC)
  out <- c(out, SE=list(rese))
} else {
  out <- c(out, BatchSE=list(coda::batchSE(out$resultMCMC, batchSize=batchSize)))
}

class(out) <- "mcmcComposite"

out <- c(out, TimeSeriesSE=list(summary(out)$statistics[,4]))

class(out) <- "mcmcComposite"

return(out)
}
