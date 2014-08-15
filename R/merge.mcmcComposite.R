#' merge.mcmcComposite Merge two mcmcComposite results
#' @title Merge two mcmcComposite results
#' @author Marc Girondot
#' @return A mcmcComposite result
#' @param x A mcmcComposite result
#' @param y A mcmcComposite result
#' @param ... Other parameters for merge function (not used)
#' @description Merge two mcmcComposite results
#' @examples
#' \dontrun{
#' data(result_mcmc_4p)
#' # Of course you must do it with two different mcmc objects !
#' result_mcmc <- merge.mcmcComposite(result_mcmc_4p, result_mcmc_4p)
#' }
#' @export

merge.mcmcComposite <- function(x, y, ...) {
  
  mcmcComposite1 <- x
  mcmcComposite2 <- y
  
  if (mcmcComposite1$parametersMCMC$n.chains != mcmcComposite2$parametersMCMC$n.chains) {
    print("MCMC results must have the same number of chains")
    return()
  }
  
  mcmcComposite <- mcmcComposite1
  for (chain in 1:mcmcComposite1$parametersMCMC$n.chains) {
  mcmcComposite$resultMCMC[[chain]] <- mcmc(rbind(mcmcComposite1$resultMCMC[[chain]], 
                                             mcmcComposite2$resultMCMC[[chain]]), 1, 
                                            mcmcComposite1$parametersMCMC$n.iter+
                                              mcmcComposite2$parametersMCMC$n.iter, 1)
  }
  
  mcmcComposite$resultLnL[[1]] <- c(mcmcComposite1$resultLnL[[1]], 
                               mcmcComposite2$resultLnL[[1]])
  
  mcmcComposite$parametersMCMC$n.iter <- mcmcComposite1$parametersMCMC$n.iter+
    mcmcComposite2$parametersMCMC$n.iter
  
  e <- as.mcmc(mcmcComposite)
  
  mcmcComposite$BatchSE <- batchSE(e)
  mcmcComposite$TimeSeriesSE <- summary(e)$statistics[,"Time-series SE"]
  return(mcmcComposite)
}