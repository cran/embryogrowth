#' STRN_MHmcmc runs the Metropolis-Hastings algorithm for STRN (Bayesian MCMC)
#' @title Metropolis-Hastings algorithm for Sexualisation Thermal Reaction Norm
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
#' @param result An object obtained after a STRN fit
#' @param parametersMCMC A set of parameters used as initial point for searching with information on priors
#' @param n.iter Number of iterations for each step
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param trace TRUE or FALSE or period, shows progress
#' @param traceML TRUE or FALSE to show ML
#' @param batchSize Number of observations to include in each batch fo SE estimation
#' @param filename If intermediate is not NULL, save intermediate result in this file
#' @param intermediate Period for saving intermediate result, NULL for no save
#' @param adaptive Should an adaptive process for SDProp be used
#' @param adaptive.lag  Lag to analyze the SDProp value in an adaptive content
#' @param adaptive.fun Function used to change the SDProp
#' @param parallel Should parallel computing for info.nests() be used
#' @param previous Previous result to be continued. Can be the filename in which intermediate results are saved.
#' @param fill Parameters to be sent to STRN().
#' @param WAIC Prepare the output for loo() and waic().
#' @description Run the Metropolis-Hastings algorithm for Sexualisation Thermal Reaction Norm.\cr
#' The number of iterations is n.iter+n.adapt+1 because the initial likelihood is also displayed.\cr
#' I recommend that thin=1 because the method to estimate SE uses resampling.\cr
#' If initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
#' To get the SE of the point estimates from \code{result_mcmc <- STRN_MHmcmc(result=try)}, use:\cr
#' \code{result_mcmc$SD}\cr
#' coda package is necessary for this function.\cr
#' The parameters intermediate and filename are used to save intermediate results every 'intermediate' iterations (for example 1000). Results are saved in a file of name filename.\cr
#' The parameter previous is used to indicate the list that has been save using the parameters intermediate and filename. It permits to continue a mcmc search.\cr
#' These options are used to prevent the consequences of computer crash or if the run is very very long and processes at time limited.\cr
#' If fill is NA, it will use the stored fill value in result.\cr
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' MedIncubation_Cc <- subset(DatabaseTSD, Species=="Caretta caretta" & 
#' RMU=="Mediterranean" & Sexed!=0)
#' Med_Cc <- tsd(males=MedIncubation_Cc$Males, 
#'              females=MedIncubation_Cc$Females, 
#'              temperatures=MedIncubation_Cc$Incubation.temperature, 
#'              par=c(P=29.5, S=-0.1))
#' plot(Med_Cc, xlim=c(25, 35))
#' males <- c(7, 0, 0, 0, 0, 5, 6, 3, 5, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' names(males) <- rev(rev(names(resultNest_4p_SSM$data))[-(1:2)])
#' sexed <- rep(10, length(males))
#' names(sexed) <- rev(rev(names(resultNest_4p_SSM$data))[-(1:2)])
#' Initial_STRN <- resultNest_4p_SSM$par[c("DHA", "DHH", "T12H")]
#' Initial_STRN <- structure(c(582.567096666926, 2194.0806711639, 3475.28414940385), 
#'                           .Names = c("DHA", "DHH", "T12H"))
#' fp <- c(Rho25=100)
#' fitSTRN <- STRN(Initial_STRN=Initial_STRN, 
#'                 EmbryoGrowthTRN=resultNest_4p_SSM, 
#'                 tsd=Med_Cc, 
#'                 embryo.stages="Caretta caretta.SCL", 
#'                 Sexed=sexed, Males=males, 
#'                 fixed.parameters=fp, 
#'                 sexratio="TSP.GrowthWeighted.STRNWeighted.sexratio")
#'                 
#' pMCMC <- TRN_MHmcmc_p(fitSTRN, accept=TRUE)
#' pMCMC[, "Density"] <- "dunif"
#' pMCMC[, "Prior2"] <- pMCMC[, "Max"]<- 10000
#' pMCMC[, "Prior1"] <- pMCMC[, "Min"] <- 1
#' outMCMC <- STRN_MHmcmc(result = fitSTRN, n.iter = 10000, parametersMCMC = pMCMC,
#'                 n.chains = 1, n.adapt = 0, thin = 1, trace = TRUE,
#'                 adaptive = TRUE, adaptive.lag = 500, 
#'                 intermediate = 1000,
#'                 filename = "intermediate_mcmcSTRN.Rdata")
#' plot(outMCMC, parameters=1)
#' plot(outMCMC, parameters=2)
#' plot(outMCMC, parameters=3)
#' 1-rejectionRate(as.mcmc(x = outMCMC))
#' 
#' }
#' @export

STRN_MHmcmc <- function(result=NULL, n.iter=10000, 
                        parametersMCMC=NULL, n.chains = 1, n.adapt = 0, 
                        thin=1, trace=NULL, traceML=FALSE, batchSize=sqrt(n.iter), 
                        adaptive=FALSE, adaptive.lag=500, 
                        adaptive.fun=function(x) {ifelse(x>0.234, 1.3, 0.7)},
                        parallel=TRUE, WAIC=TRUE, 
                        intermediate=NULL, filename="intermediate.Rdata", 
                        previous=NULL, fill = NA) {
  
  # result=NULL; n.iter=10000; parametersMCMC=NULL; n.chains = 1; n.adapt = 0; thin=1; trace=NULL; batchSize=sqrt(n.iter); adaptive=FALSE; adaptive.lag=500; adaptive.fun=function(x) {ifelse(x>0.234, 1.3, 0.7)}; intermediate=NULL; filename="intermediate.Rdata"; previous=NULL
  
  if (!is.null(fill)) {
    if (is.na(fill)) fill <- result$fill
  }
  
  if (is.character(previous)) {
    itr <- NULL
    load(previous)
    previous <- itr
    rm(itr)
    print("Continue previous mcmc run")
  } else {
    print(parametersMCMC)
  }
  
  if (!is.null(previous)) if (!is.null(trace)) previous$trace <- trace
  
  if (is.null(trace)) trace <- FALSE
  
  # 29/1/2014; Ajout de result$weight
  out <- MHalgoGen(n.iter=n.iter, 
                   parameters=parametersMCMC, 
                   parameters_name="par", 
                   n.chains = n.chains, n.adapt = n.adapt, thin=thin, 
                   trace=trace, traceML=traceML, 
                   likelihood=getFromNamespace(".STRN_fit", ns="embryogrowth"),
                   fixed.parameters=result$fixed.parameters, 
                   equation=result$equation, 
                   tsd=result$tsd, 
                   Sexed=result$Sexed, 
                   Males=result$Males, 
                   sexratio=result$sexratio, 
                   zero = 1E-9                                               , 
                   NestsResult = result$EmbryoGrowthTRN                      , 
                   embryo.stages = result$embryo.stages                      ,
                   TSP.borders = result$TSP.borders                          ,
                   TSP.begin = result$TSP.begin                              ,
                   TSP.end = result$TSP.end                                  ,
                   fill = fill                                               , 
                   WAIC = WAIC                                               ,
                   WAIC.out = WAIC                                           ,
                   n.datapoints = length(result$Sexed)                       ,
                   parallel=parallel,
                   adaptive=adaptive, 
                   adaptive.lag=adaptive.lag, 
                   adaptive.fun=adaptive.fun,
                   intermediate=intermediate, 
                   filename=filename, 
                   previous=previous)
  
  fin <- try(summary(out), silent=TRUE)
  
  if (batchSize>=n.iter/2) {
    print("batchSize cannot be larger than half the number of iterations.")
    rese <- rep(NA, dim(parametersMCMC)[1])
    names(rese) <- rownames(parametersMCMC)
    out <- c(out, SE=list(rese))
  } else {
    out <- c(out, BatchSE=list(coda::batchSE(out$resultMCMC, batchSize=batchSize)))
  }
  
  if (inherits(fin, "try-error")) { 
    lp <- rep(NA, nrow(out$parametersMCMC$parameters))
    names(lp) <- rownames(out$parametersMCMC$parameters)
    out <- c(out, TimeSeriesSE=list(lp))
    out <- c(out, SD=list(lp))
  } else {
    out <- c(out, TimeSeriesSE=list(fin$statistics[,4]))
    out <- c(out, SD=list(fin$statistics[,"SD"]))
  }
  
  out <- addS3Class(out, "mcmcComposite")
  
  return(out)
}
