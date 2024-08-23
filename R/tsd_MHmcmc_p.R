#' tsd_MHmcmc_p generates set of parameters to be used with tsd_MHmcmc()
#' @title Generates set of parameters to be used with tsd_MHmcmc()
#' @author Marc Girondot
#' @return A matrix with the parameters
#' @param result An object obtained after a tsd fit
#' @param default The default distribution for priors; can be dnorm only at that time
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
#' 
#' }
#' @export


tsd_MHmcmc_p<-function(result=stop("An output from tsd() must be provided"), 
                       default="dnorm"                                     , 
                       accept=FALSE                                        ) {
  
  default <- tolower(default)
  
  default <- match.arg(default, 
                       choices = c("dnorm", "dunif"), 
                       several.ok = FALSE)
  
  # d'abord je sors les paramtres  utiliser
  par <- result$par
  
  if (default == "dnorm") {
    
    # "P"
    # C'est seulement vrai en temperatures
    if (result$type == "temperature") {
      P <- list("dnorm", par["P"], 2, 2, 25, 35, ifelse(is.na(par["P"]), 29.5, par["P"]))
    } else {
      P <- list("dnorm", par["P"], 10, 5, 40, 100, ifelse(is.na(par["P"]), 55, par["P"]))
    }
    S <- list("dnorm", par["S"], 1, 0.5, -2, 2, ifelse(is.na(par["S"]), 0.01, par["S"]))
    K <- list("dnorm", par["K"], 3, 0.5, -20, 20, ifelse(is.na(par["K"]), 0, par["K"]))
    K1 <- list("dnorm", par["K1"], 20, 0.5, min(-100, par["K1"]-100), max(100, par["K1"]+100), ifelse(is.na(par["K1"]), 1, par["K1"]))
    K2 <- list("dnorm", par["K2"], 20, 0.5, min(-100, par["K2"]-100), max(100, par["K2"]+100), ifelse(is.na(par["K2"]), 1, par["K2"]))
    if (result$type == "temperature") {
      P_low <- list("dnorm", par["P_low"], 2, 2, 25, 35, ifelse(is.na(par["P_low"]), 29.5, par["P_low"]))
    } else {
      P_low <- list("dnorm", par["P_low"], 10, 5, 40, 100, ifelse(is.na(par["P_low"]), 55, par["P_low"]))
    }
    S_low <- list("dnorm", par["S_low"], 1, 0.5, -2, 2, ifelse(is.na(par["S_low"]), 0.01, par["S_low"]))
    K1_low <- list("dnorm", par["K1_low"], 20, 0.5, min(-100, par["K1_low"]-100), max(100, par["K1_low"]+100), ifelse(is.na(par["K1_low"]), 1, par["K1_low"]))
    K2_low <- list("dnorm", par["K2_low"], 20, 0.5, min(-100, par["K2_low"]-100), max(100, par["K2_low"]+100), ifelse(is.na(par["K2_low"]), 1, par["K2_low"]))
    if (result$type == "temperature") {
      P_high <- list("dnorm", par["P_high"], 2, 2, 25, 35, ifelse(is.na(par["P_high"]), 29.5, par["P_high"]))
    } else {
      P_high <- list("dnorm", par["P_high"], 10, 5, 40, 100, ifelse(is.na(par["P_high"]), 55, par["P_high"]))
    }
    S_high <- list("dnorm", par["S_high"], 1, 0.5, -2, 2, ifelse(is.na(par["S_high"]), 0.01, par["S_high"]))
    K1_high <- list("dnorm", par["K1_high"], 20, 0.5, min(-100, par["K1_high"]-100), max(100, par["K1_high"]+100), ifelse(is.na(par["K1_high"]), 1, par["K1_high"]))
    K2_high <- list("dnorm", par["K2_high"], 20, 0.5, min(-100, par["K2_high"]-100), max(100, par["K2_high"]+100), ifelse(is.na(par["K2_high"]), 1, par["K2_high"]))
    
    SL <- list("dnorm", par["SL"], 1, 0.2, min(-5, par["SL"] - 5), max(5, par["SL"] + 5), ifelse(is.na(par["SL"]), -1, par["SL"]))
    SH <- list("dnorm", par["SH"], 1, 0.2, min(-5, par["SH"] - 5), max(5, par["SH"] + 5), ifelse(is.na(par["SH"]), -1, par["SH"]))
    TransitionS <- list("dnorm", par["TransitionS"], 0, 0.4, -5, 5, ifelse(is.na(par["TransitionS"]), 0, par["TransitionS"]))
    
    SL_low <- list("dnorm", par["SL_low"], 1, 0.2, min(-5, par["SL_low"] - 5), max(5, par["SL_low"] + 5), ifelse(is.na(par["SL_low"]), -1, par["SL_low"]))
    SH_low <- list("dnorm", par["SH_low"], 1, 0.2, min(-5, par["SH_low"] - 5), max(5, par["SH_low"] + 5), ifelse(is.na(par["SH_low"]), -1, par["SH_low"]))
    TransitionS_low <- list("dnorm", par["TransitionS_low"], 0, 0.4, -5, 5, ifelse(is.na(par["TransitionS_low"]), 0, par["TransitionS_low"]))
    SL_high <- list("dnorm", par["SL_high"], 1, 0.2, min(-5, par["SL_high"] - 5), max(5, par["SL_high"] + 5), ifelse(is.na(par["SL_high"]), -1, par["SL_high"]))
    SH_high <- list("dnorm", par["SH_high"], 1, 0.2, min(-5, par["SH_high"] - 5), max(5, par["SH_high"] + 5), ifelse(is.na(par["SH_high"]), -1, par["SH_high"]))
    TransitionS_high <- list("dnorm", par["TransitionS_high"], 0, 0.4, -5, 5, ifelse(is.na(par["TransitionS_high"]), 0, par["TransitionS_high"]))
  } else {
    if (result$type == "temperature") {
      P <- list("dunif", 25, 35, 2, 25, 35, ifelse(is.na(par["P"]), 29.5, par["P"]))
    } else {
      P <- list("dunif", 40, 100, 2, 40, 100, ifelse(is.na(par["P"]), 55, par["P"]))
    }
    S <- list("dunif", -2, 2, 0.2, -2, 2, ifelse(is.na(par["S"]), 0.01, par["S"]))
    K <- list("dunif", -20, 20, 0.2, -20, 20, ifelse(is.na(par["K"]), 0, par["K"]))
    K1 <- list("dunif", min(-100, par["K1"]-100), max(100, par["K1"]+100), 0.2, min(-100, par["K1"]-100), max(100, par["K1"]+100), ifelse(is.na(par["K1"]), 1, par["K1"]))
    K2 <- list("dunif", min(-100, par["K2"]-100), max(100, par["K2"]+100), 0.2, min(-100, par["K2"]-100), max(100, par["K2"]+100), ifelse(is.na(par["K2"]), 1, par["K2"]))
    if (result$type == "temperature") {
      P_low <- list("dunif", 25, 35, 2, 25, 35, ifelse(is.na(par["P_low"]), 29.5, par["P_low"]))
    } else {
      P_low <- list("dunif", 40, 100, 2, 40, 100, ifelse(is.na(par["P_low"]), 55, par["P_low"]))
    }
    S_low <- list("dunif", -2, 2, 0.2, -2, 2, ifelse(is.na(par["S_low"]), 0.01, par["S_low"]))
    K1_low <- list("dunif", min(-100, par["K1_low"]-100), max(100, par["K1_low"]+100), 0.2, min(-100, par["K1_low"]-100), max(100, par["K1_low"]+100), ifelse(is.na(par["K1_low"]), 1, par["K1_low"]))
    K2_low <- list("dunif", min(-100, par["K2_low"]-100), max(100, par["K2_low"]+100), 0.2, min(-100, par["K2_low"]-100), max(100, par["K2_low"]+100), ifelse(is.na(par["K2_low"]), 1, par["K2_low"]))
    if (result$type == "temperature") {
      P_high <- list("dunif", 25, 35, 2, 25, 35, ifelse(is.na(par["P_high"]), 29.5, par["P_high"]))
    } else {
      P_high <- list("dunif", 40, 100, 2, 40, 100, ifelse(is.na(par["P_high"]), 55, par["P_high"]))
    }
    S_high <- list("dunif", -2, 2, 0.2, -2, 2, ifelse(is.na(par["S_high"]), 0.01, par["S_high"]))
    K1_high <- list("dunif", min(-100, par["K1_high"]-100), max(100, par["K1_high"]+100), 0.2, min(-100, par["K1_high"]-100), max(100, par["K1_high"]+100), ifelse(is.na(par["K1_high"]), 1, par["K1_high"]))
    K2_high <- list("dunif", min(-100, par["K2_high"]-100), max(100, par["K2_high"]+100), 0.2, min(-100, par["K2_high"]-100), max(100, par["K2_high"]+100), ifelse(is.na(par["K2_high"]), 1, par["K2_high"]))
    
    SL <- list("dunif", min(-5, par["SL"] - 5), max(5, par["SL"] + 5), 0.2, min(-5, par["SL"] - 5), max(5, par["SL"] + 5), ifelse(is.na(par["SL"]), -1, par["SL"]))
    SH <- list("dunif", min(-5, par["SH"] - 5), max(5, par["SH"] + 5), 0.2, min(-5, par["SH"] - 5), max(5, par["SH"] + 5), ifelse(is.na(par["SH"]), -1, par["SH"]))
    TransitionS <- list("dunif", -5, 5, 0.2, -5, 5, ifelse(is.na(par["TransitionS"]), 0, par["TransitionS"]))
    
    SL_low <- list("dunif", min(-5, par["SL_low"] - 5), max(5, par["SL_low"] + 5), 0.2, min(-5, par["SL_low"] - 5), max(5, par["SL_low"] + 5), ifelse(is.na(par["SL_low"]), -1, par["SL_low"]))
    SH_low <- list("dunif", min(-5, par["SH_low"] - 5), max(5, par["SH_low"] + 5), 0.2, min(-5, par["SH_low"] - 5), max(5, par["SH_low"] + 5), ifelse(is.na(par["SH_low"]), -1, par["SH_low"]))
    TransitionS_low <- list("dunif", -5, 5, 0.2, -5, 5, ifelse(is.na(par["TransitionS_low"]), 0, par["TransitionS_low"]))
    SL_high <- list("dunif", min(-5, par["SL_high"] - 5), max(5, par["SL_high"] + 5), 0.2, min(-5, par["SL_high"] - 5), max(5, par["SL_high"] + 5), ifelse(is.na(par["SL_high"]), -1, par["SL_high"]))
    SH_high <- list("dunif", min(-5, par["SH_high"] - 5), max(5, par["SH_high"] + 5), 0.2, min(-5, par["SH_high"] - 5), max(5, par["SH_high"] + 5), ifelse(is.na(par["SH_high"]), -1, par["SH_high"]))
    TransitionS_high <- list("dunif", -5, 5, -0.2, 5, 5, ifelse(is.na(par["TransitionS_high"]), 0, par["TransitionS_high"]))
  }
  parameters <- rbind(P=P, S=S, K=K, K1=K1, K2=K2, 
                 SL=SL, SH=SH, TransitionS=TransitionS, 
                 SL_low=SL_low, SH_low=SH_low, TransitionS_low=TransitionS_low, 
                 SL_high= SL_high, SH_high= SH_high, TransitionS_high=TransitionS_high, 
                 P_low=P_low, S_low=S_low, K1_low=K1_low, K2_low=K2_low, 
                 P_high=P_high, S_high=S_high, K1_high=K1_high, K2_high=K2_high)
  
  # prencours <- NULL
  # 
  # for (i in 1:length(par)) {
  #   prencours <- c(prencours, priors[[names(par)[i]]])
  # }
  # 
  # parametersMCMC <- matrix(prencours, ncol=7, byrow=T)
  colnames(parameters) <- c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", "Init")
  # rownames(parametersMCMC)<-names(par)
  
  # parametersMCMC <- as.data.frame(parametersMCMC, stringsAsFactors = FALSE)
  # 
  # for (i in c("Prior1", "Prior2", "SDProp", "Min", "Max", "Init"))
  #   parametersMCMC[,i] <- as.numeric(parametersMCMC[,i])
  # 
  # parameters <- parametersMCMC
  
  if (accept) {
    
    return(parameters[names(par), ])
    
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
