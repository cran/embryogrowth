#' HatchingSuccess.fit fits a hatching success model to data
#' @title Fit a hatching success model to data using maximum likelihood
#' @author Marc Girondot
#' @return Return a object of class HatchingSuccess
#' @param par A set of parameters.
#' @param data A dataset in a data.frame with a least three columns: Incubation.temperature, Hatched and NotHatched
#' @param fixed.parameters A set of parameters that must not be fitted.
#' @param column.Incubation.temperature Name of the column with incubation temperatures
#' @param column.Hatched Name of the column with hatched number
#' @param column.NotHatched Name of the column with not hatched number
#' @param hessian Should Hessian matrix be estimated?
#' @description Set of functions to study the hatching success.\cr
#' The first version of the model was published in:\cr
#' Laloë, J.-O., Monsinjon, J., Gaspar, C., Touron, M., Genet, Q., Stubbs, J., Girondot, M. 
#' & Hays, G.C. (2020) Production of male hatchlings at a remote South Pacific green sea turtle 
#' rookery: conservation implications in a female-dominated world. Marine Biology, 167, 70.\cr
#' The version available here is enhanced by using a double flexit model rather than a double 
#' logistic model. The flexit model is described here:\cr
#' Abreu-Grobois, F.A., Morales-Mérida, B.A., Hart, C.E., Guillon, J.-M., Godfrey, M.H., 
#' Navarro, E. & Girondot, M. (2020) Recent advances on the estimation of the thermal 
#' reaction norm for sex ratios. PeerJ, 8, e8451.\cr
#' @family Hatching success
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' totalIncubation_Cc <- subset(DatabaseTSD, 
#'                              Species=="Caretta caretta" & 
#'                                Note != "Sinusoidal pattern" & 
#'                                !is.na(Total) & Total != 0)
#' 
#' par <- c(S.low=0.5, S.high=0.3, 
#'          P.low=25, deltaP=10, MaxHS=0.8)
#'          
#' g.logistic <- HatchingSuccess.fit(par=par, data=totalIncubation_Cc)
#'          
#' HatchingSuccess.lnL(par=g.logistic$par, data=totalIncubation_Cc)
#' 
#' plot(g.logistic)
#' 
#' par <- c(S.low=0.5, S.high=0.3, 
#'          P.low=25, deltaP=10, 
#'          K1.low=1, K2.low=-1, K1.high=1, K2.high=-1, 
#'          MaxHS=0.8)
#' 
#' g.flexit <- HatchingSuccess.fit(par=par, data=totalIncubation_Cc)
#' 
#' HatchingSuccess.lnL(par=g.flexit$par, data=totalIncubation_Cc)
#' 
#' compare_AICc(logistic=g.logistic, flexit=g.flexit)
#' plot(x=g.logistic, what = c("observations", "ML", "CI"), replicates=10000)
#' 
#' pMCMC <- HatchingSuccess.MHmcmc_p(result = g.logistic, accept = TRUE)
#' MCMC <- HatchingSuccess.MHmcmc(result = g.logistic, parametersMCMC = pMCMC,
#'                             n.iter = 100000, 
#'                            adaptive = TRUE)
#' 
#' plot(MCMC, parameters = "S.low")
#' plot(MCMC, parameters = "S.high")
#' plot(MCMC, parameters = "P.low")
#' plot(MCMC, parameters = "deltaP")
#' plot(MCMC, parameters = "MaxHS")
#' 
#' plot(x=g.logistic, what = c("observations", "ML", "CI"), 
#'         replicates=10000, resultmcmc=MCMC)
#' 
#' ####### Exemple with Chelonia mydas
#' totalIncubation_Cm <- subset(DatabaseTSD, 
#' Species=="Chelonia mydas" & 
#'   Note != "Sinusoidal pattern" & 
#'   !is.na(Total) & Total != 0 & 
#'   !is.na(NotHatched) & !is.na(Hatched))
#'   
#'   totalIncubation_Cm$NotHatched <- totalIncubation_Cm$NotHatched + 
#'   ifelse(!is.na(totalIncubation_Cm$Undeveloped), totalIncubation_Cm$Undeveloped, 0)
#'   
#'   
#'   plot(x=totalIncubation_Cm$Incubation.temperature, 
#'   y=totalIncubation_Cm$Hatched/totalIncubation_Cm$Total, bty="n", las=1, 
#'   xlab="Constant incubation temperature", ylab="Proportion of hatching")
#' 
#' par <- c(S.low=0.5, S.high=0.3, 
#'          P.low=25, deltaP=10, MaxHS=0.8)
#' 
#' g.logistic <- HatchingSuccess.fit(par=par, data=totalIncubation_Cm)
#' plot(g.logistic)
#' 
#' pMCMC <- HatchingSuccess.MHmcmc_p(g.logistic, accept=TRUE)
#' mcmc <- HatchingSuccess.MHmcmc(result=g.logistic, parameters = pMCMC, 
#'                   adaptive=TRUE, n.iter=100000, trace=1000)
#' par <- as.parameters(mcmc)
#' par <- as.parameters(mcmc, index="median")
#' 
#' plot(mcmc, parameters=c("P.low"))
#' plot(mcmc, parameters=c("deltaP"))
#' plot(mcmc, parameters=c("S.low"))
#' plot(mcmc, parameters=c("S.high"))
#' plot(mcmc, parameters=c("MaxHS"))
#' 
#' plot(g.logistic, resultmcmc=mcmc, what = c("observations", "CI"))
#' }
#' @export


HatchingSuccess.fit <- function(par=NULL, data=stop("data must be provided"), 
                                fixed.parameters=NULL, 
                                column.Incubation.temperature="Incubation.temperature",
                                column.Hatched="Hatched",
                                column.NotHatched="NotHatched", 
                                hessian=TRUE) {
  
  # par=NULL; data=NULL; fixed.parameters=NULL; column.Incubation.temperature="Incubation.temperature"; column.Hatched="Hatched"; column.NotHatched="NotHatched"; hessian=FALSE
  
  if (is.null(par)) par <- c('S.low' = 1.426059810893495, 
                             'S.high' = 0.10953915294233023, 
                             'P.low' = 15.9139408256423156, 
                             'P.high' = 32.494793698343166, 
                             'MaxHS' = 0.99)
  
  lw <- c('S.low' = 0, 'S.high' = 0, 'P.low' = 10, 'P.high' = 20, 'MaxHS' = 0, K1.low=-10, K1.high=-10, K2.low=-10, K2.high=-10)[names(par)]
  hg <- c('S.low' = 20, 'S.high' = 20, 'P.low' = 50, 'P.high' = 100, 'MaxHS' = 1, K1.low=+10, K1.high=+10, K2.low=+10, K2.high=+10)[names(par)]
  
  g <- suppressWarnings(optim(par=par, fn=HatchingSuccess.lnL, 
                              fixed.parameters=fixed.parameters, 
                              column.Incubation.temperature=column.Incubation.temperature,
                              column.Hatched=column.Hatched,
                              column.NotHatched=column.NotHatched,
                              method="L-BFGS-B", 
                              lower=lw, 
                              upper=hg, 
                              data=data, hessian = FALSE))
  
  N <- sum(data[, column.Hatched]) + sum(data[, column.NotHatched])
  
  # g$par["S.low"] <- abs(g$par["S.low"])
  # g$par["S.high"] <- abs(g$par["S.high"])
  
  g <- suppressWarnings(optim(par=g$par, fn=HatchingSuccess.lnL, 
                              fixed.parameters=fixed.parameters, 
                              column.Incubation.temperature=column.Incubation.temperature,
                              column.Hatched=column.Hatched,
                              column.NotHatched=column.NotHatched,
                              method="L-BFGS-B", 
                              lower=lw, 
                              upper=hg, 
                              data=data, hessian = hessian))
  
  g$fixed.parameters <- fixed.parameters
  if (hessian) {
    g$SE <- SEfromHessian(g$hessian)
  }
  g$AIC <- 2*g$value+2*length(g$par)
  # Correction le 23/11/2020
  g$AICc <- g$AIC +(2*length(par)*(length(par)+1))/(sum(N)-length(par)-1)
  # Correction le 23/11/2020
  g$BIC <- 2*g$value+ length(par)*log(sum(N))
  
  g$data <- data
  g$column.Incubation.temperature <- column.Incubation.temperature
  g$column.Hatched <- column.Hatched
  g$column.NotHatched <- column.NotHatched
  
  g <- addS3Class(g, "HatchingSuccess")
  return(g)
}


