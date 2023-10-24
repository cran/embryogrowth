#' tsd_MHmcmc_p generates set of parameters to be used with tsd_MHmcmc()
#' @title Generates set of parameters to be used with tsd_MHmcmc()
#' @author Marc Girondot
#' @return A matrix with the parameters
#' @param result An object obtained after a tsd fit
#' @param accept If TRUE, the script does not wait user information
#' @description Interactive script used to generate set of parameters to be 
#' used with tsd_MHmcmc().
#' @family Functions for temperature-dependent sex determination
#' @examples 
#' \dontrun{
#' library(embryogrowth)
#' eo <- subset(DatabaseTSD, Species=="Emys orbicularis", c("Males", "Females", 
#'                                        "Incubation.temperature"))
#' eo_logistic <- with(eo, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature))
#' pMCMC <- tsd_MHmcmc_p(eo_logistic, accept=TRUE)
#' # Take care, it can be very long
#' result_mcmc_tsd <- tsd_MHmcmc(result=eo_logistic, 
#' 		parametersMCMC=pMCMC, n.iter=10000, n.chains = 1,  
#' 		n.adapt = 0, thin=1, trace=FALSE, adaptive=TRUE)
#' # summary() permits to get rapidly the standard errors for parameters
#' summary(result_mcmc_tsd)
#' plot(result_mcmc_tsd, parameters="S", scale.prior=TRUE, xlim=c(-3, 3), las=1)
#' plot(result_mcmc_tsd, parameters="P", scale.prior=TRUE, xlim=c(25, 35), las=1)
#' 
#' eo_flexit <- with(eo, tsd(males=Males, females=Females,
#'                                  parameters.initial=c(eo_logistic$par["P"], 
#'                                  1/(4*eo_logistic$par["S"]), 
#'                                  K1=1, K2=1), 
#'                                  temperatures=Incubation.temperature,
#'                                  equation="flexit", replicate.CI=NULL))
#' pMCMC <- tsd_MHmcmc_p(eo_flexit, accept=TRUE)
#' result_mcmc_tsd <- tsd_MHmcmc(result=eo_flexit, 
#' 		parametersMCMC=pMCMC, n.iter=10000, n.chains = 1,  
#' 		n.adapt = 0, thin=1, trace=FALSE, adaptive=TRUE)
#' # summary() permits to get rapidly the standard errors for parameters
#' summary(result_mcmc_tsd)
#' plot(result_mcmc_tsd, parameters="S", scale.prior=TRUE, xlim=c(-3, 3), las=1)
#' plot(result_mcmc_tsd, parameters="P", scale.prior=TRUE, xlim=c(25, 35), las=1)
#' plot(result_mcmc_tsd, parameters="K1", scale.prior=TRUE, xlim=c(-10, 10), las=1)
#' plot(result_mcmc_tsd, parameters="K2", scale.prior=TRUE, xlim=c(-10, 10), las=1)
#' 
#' plot(eo_flexit, resultmcmc = result_mcmc_tsd)
#' }
#' @export


tsd_MHmcmc_p<-function(result=stop("An output from tsd() must be provided"), 
                       accept=FALSE) {
  
  # d'abord je sors les paramtres  utiliser
  
  par <- result$par
  
  # "P"
  # C'est seulement vrai en temperatures
  if (result$type == "temperature") {
    P <- c("dnorm", par["P"], 2, 2, 25, 35, ifelse(is.na(par["P"]), 29.5, par["P"]))
  } else {
    P <- c("dnorm", par["P"], 10, 5, 40, 100, ifelse(is.na(par["P"]), 55, par["P"]))
  }
  S <- c("dnorm", par["S"], 1, 0.5, -2, 2, ifelse(is.na(par["S"]), 0.01, par["S"]))
  K <- c("dnorm", par["K"], 3, 0.5, -20, 20, ifelse(is.na(par["K"]), 0, par["K"]))
  K1 <- c("dnorm", par["K1"], 20, 0.5, min(-100, par["K1"]-100), max(100, par["K1"]+100), ifelse(is.na(par["K1"]), 1, par["K1"]))
  K2 <- c("dnorm", par["K2"], 20, 0.5, min(-100, par["K2"]-100), max(100, par["K2"]+100), ifelse(is.na(par["K2"]), 1, par["K2"]))
  if (result$type == "temperature") {
    P_low <- c("dnorm", par["P_low"], 2, 2, 25, 35, ifelse(is.na(par["P_low"]), 29.5, par["P_low"]))
  } else {
    P_low <- c("dnorm", par["P_low"], 10, 5, 40, 100, ifelse(is.na(par["P_low"]), 55, par["P_low"]))
  }
  S_low <- c("dnorm", par["S_low"], 1, 0.5, -2, 2, ifelse(is.na(par["S_low"]), 0.01, par["S_low"]))
  K1_low <- c("dnorm", par["K1_low"], 20, 0.5, min(-100, par["K1_low"]-100), max(100, par["K1_low"]+100), ifelse(is.na(par["K1_low"]), 1, par["K1_low"]))
  K2_low <- c("dnorm", par["K2_low"], 20, 0.5, min(-100, par["K2_low"]-100), max(100, par["K2_low"]+100), ifelse(is.na(par["K2_low"]), 1, par["K2_low"]))
  if (result$type == "temperature") {
    P_high <- c("dnorm", par["P_high"], 2, 2, 25, 35, ifelse(is.na(par["P_high"]), 29.5, par["P_high"]))
  } else {
    P_high <- c("dnorm", par["P_high"], 10, 5, 40, 100, ifelse(is.na(par["P_high"]), 55, par["P_high"]))
  }
  S_high <- c("dnorm", par["S_high"], 1, 0.5, -2, 2, ifelse(is.na(par["S_high"]), 0.01, par["S_high"]))
  K1_high <- c("dnorm", par["K1_high"], 20, 0.5, min(-100, par["K1_high"]-100), max(100, par["K1_high"]+100), ifelse(is.na(par["K1_high"]), 1, par["K1_high"]))
  K2_high <- c("dnorm", par["K2_high"], 20, 0.5, min(-100, par["K2_high"]-100), max(100, par["K2_high"]+100), ifelse(is.na(par["K2_high"]), 1, par["K2_high"]))
  
  
  priors <- list(P=P, S=S, K=K, K1=K1, K2=K2, 
                 P_low=P_low, S_low=S_low, K1_low=K1_low, K2_low=K2_low, 
                 P_high=P_high, S_high=S_high, K1_high=K1_high, K2_high=K2_high)
  
  prencours <- NULL
  
  for (i in 1:length(par)) {
    prencours <- c(prencours, priors[[names(par)[i]]])
  }
  
  parametersMCMC <- matrix(prencours, ncol=7, byrow=T)
  colnames(parametersMCMC) <- c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", "Init")
  rownames(parametersMCMC)<-names(par)
  
  parametersMCMC <- as.data.frame(parametersMCMC, stringsAsFactors = FALSE)
  
  for (i in c("Prior1", "Prior2", "SDProp", "Min", "Max", "Init"))
    parametersMCMC[,i] <- as.numeric(parametersMCMC[,i])
  
  parameters <- parametersMCMC
  
  if (accept) {
    
    return(parameters)
    
  } else {
    
    repeat {
      
      cat("Proposition:\n")
      print(parameters)
      cat("Name of the parameter to change or Enter to quit:\n")
      f<-scan(nmax=1, quiet=TRUE, what=character())
      
      if (length(f)==0) f <- "q"
      
      if (f=="q") {
        return(parameters)
        
      } else {
        
        variable <- which(f==names(par))
        if (length(variable)==0) {
          cat("The parameter does not exist:\n")
        } else {
          print(variable)
          cat(paste("Change for the parameter ",names(par)[variable],":\n",sep=""))
          
          cat(paste("Distribution of the prior (Enter for default ",parameters[variable, "Density"], "):", sep=""))
          density<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(density)!=0) { parameters[variable, "Density"] <- density } else { density <- parameters[variable, "Density"] }
          
          if (density == "dunif") {
            
            cat(paste("Distribution of the prior, Minimum (Enter for default ", parameters[variable, "Prior1"], "):", sep=""))
            f<-scan(nmax=1, quiet=TRUE, what=character())
            if (length(f)!=0) parameters[variable, "Prior1"] <- f
            cat(paste("Distribution of the prior, Maximum (Enter for default ", parameters[variable, "Prior2"], "):", sep=""))
            f<-scan(nmax=1, quiet=TRUE, what=character())
            if (length(f)!=0) parameters[variable, "Prior2"] <- f
            
          } else {
            
            if (density == "dnorm") {
              
              cat(paste("Distribution of the prior, Mean (Enter for default ", parameters[variable, "Prior1"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior1"] <- f
              cat(paste("Distribution of the prior, Standard deviation (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior2"] <- f
              
            } else {
              
              cat(paste("Distribution of the prior, value 1 (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior1"] <- f
              cat(paste("Distribution of the prior, value 2 (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
              f<-scan(nmax=1, quiet=TRUE, what=character())
              if (length(f)!=0) parameters[variable, "Prior2"] <- f
              
            }
          }
          
          
          cat(paste("SD of new proposition (Enter for default ",parameters[variable, "SDProp"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "SDProp"] <- f
          cat(paste("Minimum for the parameter (default ",parameters[variable, "Min"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Min"] <- f
          cat(paste("Maximum for the parameter (Enter for default ",parameters[variable, "Max"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Max"] <- f
          cat(paste("Initial value (Enter for default ",parameters[variable, "Init"], "):", sep=""))
          f<-scan(nmax=1, quiet=TRUE, what=character())
          if (length(f)!=0) parameters[variable, "Init"] <- f
        }
        
      }
    }
    
  }
  
}
