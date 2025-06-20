#' GRTRN_MHmcmc runs the Metropolis-Hastings algorithm for data (Bayesian MCMC)
#' @title Metropolis-Hastings algorithm for Embryo Growth Rate Thermal Reaction Norm
#' @author Marc Girondot
#' @return A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
#' @param n.iter Number of iterations for each step
#' @param parametersMCMC A set of parameters used as initial point for searching with information on priors
#' @param result An object obtained after a SearchR fit
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param trace TRUE or FALSE or period, shows progress
#' @param traceML TRUE or FALSE to show ML
#' @param WAIC Should WAIC data been recorded?
#' @param parallel If true, try to use several cores using parallel computing
#' @param filename If intermediate is not NULL, save intermediate result in this file
#' @param intermediate Period for saving intermediate result, NULL for no save
#' @param adaptive Should an adaptive process for SDProp be used
#' @param adaptive.lag  Lag to analyze the SDProp value in an adaptive content
#' @param adaptive.fun Function used to change the SDProp
#' @param previous Previous result to be continued. Can be the filename in which intermediate results are saved.
#' @description Run the Metropolis-Hastings algorithm for data.\cr
#' The number of iterations is \code{n.iter+n.adapt+1} because the initial likelihood is also displayed.\cr
#' I recommend that thin=1 because the method to estimate SE uses resampling.\cr
#' If initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
#' To get the SE of the point estimates from \code{result_mcmc <- GRTRN_MHmcmc(result=try)}, use:\cr
#' \code{result_mcmc$SD}\cr
#' \code{coda} package is necessary for this function.\cr
#' The parameters \code{intermediate} and \code{filename} are used to save intermediate results every 'intermediate' iterations (for example 1000). Results are saved in a file named \code{filename}.\cr
#' The parameter previous is used to indicate the list that has been save using the parameters intermediate and filename. It permits to continue a mcmc search.\cr
#' These options are used to prevent the consequences of computer crash or if the run is very very long and processes with user limited time.\cr
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "T12H", "DHA",  "DHH", "DHL", "Rho25"
#' ############################################################################
#' pfixed <- c(rK=1.208968)
#' M0 = 0.3470893 
#' ############################################################################
#' # 4 parameters
#' ############################################################################
#' x c('DHA' = 109.31113503282113, 'DHH' = 617.80695919563857, 
#'     'T12H' = 306.38890489505093, 'Rho25' = 229.37265815800225)
#'     
#' resultNest_4p_SSM <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, integral=integral.Gompertz, M0=M0, 
#' 	hatchling.metric=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p_SSM)
#' plot(resultNest_4p_SSM, xlim=c(0,70), ylimT=c(22, 32), ylimS=c(0,45), series=1, 
#' embryo.stages="Caretta caretta.SCL")
#' ############################################################################
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p_SSM, accept=TRUE)
#' # Take care, it can be very long; several days
#' resultNest_mcmc_4p_SSM <- GRTRN_MHmcmc(result=resultNest_4p_SSM , 
#'                                        adaptive = TRUE          , 
#'                                        WAIC=TRUE                , 
#' 		                                    parametersMCMC=pMCMC     , 
#' 		                                    n.iter=10000            , 
#' 		                                    n.chains = 1             ,  
#' 		                                    n.adapt = 0              , 
#' 		                                    thin=1                   , 
#' 		                                    trace=TRUE               )
#' # The SDProp in pMCMC at the beginning of the Markov chain can 
#' # be considered as non-optimal. However, the values at the end 
#' # of the Markoc chain are better due to the use of the option
#' # adaptive = TRUE. Then a good strategy is to run again the MCMC
#' # with this final set:
#' pMCMC[, "SDProp"] <- resultNest_mcmc_4p_SSM$parametersMCMC$SDProp.end
#' # Also, take the set of parameters fitted as maximim likelihood
#' # as initial set of value
#' pMCMC[, "Init"] <- as.parameters(resultNest_mcmc_4p_SSM)
#' # Then I run again the MCMC; it will ensure to get the optimal distribution
#' resultNest_mcmc_4p_SSM <- GRTRN_MHmcmc(result=resultNest_4p_SSM, 
#'                                        adaptive = TRUE         ,
#' 		                                    parametersMCMC=pMCMC    , 
#' 		                                    n.iter=10               , 
#' 		                                    n.chains = 1            ,  
#' 		                                    n.adapt = 0             , 
#' 		                                    thin=1                  , 
#' 		                                    trace=TRUE              )
#' 		
#' data(resultNest_mcmc_4p_SSM)
#' out <- as.mcmc(resultNest_mcmc_4p_SSM)
#' # This out can be used with coda package
#' # Test for stationarity and length of chain
#' require(coda)
#' heidel.diag(out)
#' raftery.diag(out)
#' # plot() can use the direct output of GRTRN_MHmcmc() function.
#' plot(resultNest_mcmc_4p_SSM, parameters=1, xlim=c(0,550))
#' plot(resultNest_mcmc_4p_SSM, parameters=3, xlim=c(290,320))
#' # summary() permits to get rapidly the standard errors for parameters
#' # They are store in the result also.
#' se <- result_mcmc_4p_SSM$SD
#' # the confidence interval is better estimated by:
#' apply(out[[1]], 2, quantile, probs=c(0.025, 0.975))
#' # The use of the intermediate method is as followed;
#' # Here the total mcmc iteration is 10000, but every 1000, intermediate
#' # results are saved in file intermediate1000.Rdata:
#' resultNest_mcmc_4p_SSM <- GRTRN_MHmcmc(result=resultNest_4p_SSM, 
#' parametersMCMC=pMCMC, n.iter=10000, n.chains = 1,  
#' n.adapt = 0, thin=1, trace=TRUE, 
#' intermediate=1000, filename="intermediate1000.Rdata")
#' # If run has been stopped for any reason, it can be resumed with:
#' resultNest_mcmc_4p_SSM <- GRTRN_MHmcmc(previous="intermediate1000.Rdata")
#' # Example to use of the epsilon parameter to get confidence level
#' resultNest_4p_epsilon <- resultNest_4p
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
#' @export

GRTRN_MHmcmc <- function(result=NULL                                           , 
                         n.iter=10000                                          , 
                         parametersMCMC=NULL                                   , 
                         n.chains = 1                                          , 
                         n.adapt = 0                                           , 
                         thin=1                                                , 
                         trace=NULL                                            , 
                         traceML=FALSE                                         , 
                         WAIC=TRUE                                             , 
                         parallel=TRUE                                         , 
                         adaptive=FALSE                                        , 
                         adaptive.lag=500                                      , 
                         adaptive.fun=function(x) {ifelse(x>0.234, 1.3, 0.7)}  ,
                         intermediate=NULL                                     , 
                         filename="intermediate.Rdata"                         , 
                         previous=NULL                                         ) {
  
  # result=NULL; n.iter=10000; parametersMCMC=NULL; n.chains = 1; n.adapt = 0; thin=1; trace=FALSE; batchSize=sqrt(n.iter); parallel=TRUE; intermediate=NULL; filename="intermediate.Rdata"; previous=NULL
  # previous <- "/Users/marc/Dropbox/DropBoxPerso/_Package_HelpersMG/previous_result_mcmc_AtlMed_simplify_2.RData"
  # result=resultNest_4p; parametersMCMC=pMCMC; n.iter=100; n.chains = 1; n.adapt = 0; thin=1; trace=TRUE; intermediate = 10
  
  # result=medFit_4p;parametersMCMC=pMCMC; n.iter=1000; n.chains = 1;n.adapt = 0; thin=1; trace=TRUE
  
  if (is.character(previous)) {
    itr <- NULL
    load(previous)
    previous <- itr
    rm(itr)
  }
  
  mc.cores <- getOption("mc.cores", detectCores())
  message(paste0("I will use ", as.character(mc.cores)," cores."))
  
  
  if (is.list(previous)) {
    
    if (!is.null(trace)) previous$trace <- trace
    
    print("Continue previous mcmc run")
    # 29/1/2014; Ajout de result$weight
    out <- MHalgoGen(previous=previous)
    
  } else {
    if (is.null(trace)) trace <- FALSE
    print(parametersMCMC)
    # 29/1/2014; Ajout de result$weight
    # n.iter=n.iter; parameters=parametersMCMC; n.chains = n.chains; n.adapt = n.adapt; thin=thin; trace=trace; data=list(data=result$data, integral=result$integral, hatchling.metric=result$hatchling.metric, M0=result$M0, fixed.parameters=result$fixed.parameters, weight=result$weight); likelihood=getFromNamespace(".fonctionMCMC", ns="embryogrowth"); intermediate=intermediate; filename=filename; previous=previous
    # 3/12/2015 j'avais un data=list(data=result$data, integral=result$integral, 
    # hatchling.metric=result$hatchling.metric, M0=result$M0, fixed.parameters=result$fixed.parameters, 
    # weight=result$weight)
    out <- MHalgoGen(n.iter=n.iter                                                  , 
                     parameters=parametersMCMC                                      , 
                     n.chains = n.chains                                            , 
                     n.adapt = n.adapt                                              , 
                     thin=thin                                                      , 
                     trace=trace                                                    , 
                     traceML=traceML                                                , 
                     intermediate=intermediate                                      , 
                     filename=filename                                              , 
                     previous=previous                                              , 
                     parallel=parallel                                              , 
                     adaptive=adaptive                                              , 
                     adaptive.lag=adaptive.lag                                      , 
                     adaptive.fun=adaptive.fun                                      ,
                     temperatures=result$data                                       , 
                     integral=result$integral                                       , 
                     hatchling.metric=result$hatchling.metric                       , 
                     M0=result$M0                                                   , 
                     fixed.parameters=result$fixed.parameters                       , 
                     WAIC=WAIC                                                      ,
                     WAIC.out=WAIC                                                      ,
                     n.datapoints = unname(result$data$IndiceT["NbTS"])  ,
                     weight=result$weight                                           , 
                     out="Likelihood"                                               , 
                     progressbar=FALSE                                              , 
                     warnings=FALSE                                                 , 
                     likelihood=getFromNamespace("info.nests", ns = "embryogrowth") )
    
  }
  
  out <- addS3Class(out, "mcmcComposite")
  
  fin <- try(summary(out), silent=TRUE)
  
  if (inherits(fin, "try-error")) { 
    lp <- rep(NA, nrow(out$parametersMCMC$parameters))
    names(lp) <- rownames(out$parametersMCMC$parameters)
    out <- c(out, TimeSeriesSE=list(lp))
    out <- c(out, SD=list(lp))
  } else {
    if (is.null(nrow(fin$statistics))) {
      out <- c(out, TimeSeriesSE=list(fin$statistics[4]))
      out <- c(out, SD=list(fin$statistics["SD"]))
    } else {
      out <- c(out, TimeSeriesSE=list(fin$statistics[,4]))
      out <- c(out, SD=list(fin$statistics[,"SD"]))
    }
  }
  
  out <- addS3Class(out, "mcmcComposite")
  
  return(out)
}
