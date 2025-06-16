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
                       accept=TRUE                                        ) {
  # d'abord je sors les paramtres  utiliser
  par <- result$par
  
  default <- tolower(default)
  
  if (length(par) != length(default)) {
    default <- rep(default, length(par))[seq_along(par)]
  }
  
  default <- match.arg(default, 
                       choices = c("dnorm", "dunif"), 
                       several.ok = TRUE)
  
  minx <- ifelse(result$type == "temperature", 25, 30)
  maxx <- ifelse(result$type == "temperature", 35, 80)
  sdpx <- ifelse(result$type == "temperature", 2, 10)
  
  equation <- result$equation
  
  dpriors_norm <- NULL
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=par["P"], Prior2=sdpx, SDProp=sdpx, Min=min(c(minx, par["P"]-5)), Max=max(c(maxx, par["P"]+5)), Init=par["P"], row.names ="P"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["S"]), 0, par["S"])), Prior2=0.5, SDProp=0.5, Min=min(c(-2, ifelse(is.na(par["S"]), 0,par["S"])-2)), Max=max(c(2, ifelse(is.na(par["S"]), 0,par["S"])+2)), Init=unname(ifelse(is.na(par["S"]), 0,par["S"])), row.names = "S"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["K"]), 0, par["K"])), Prior2=3, SDProp=0.5, Min=min(c(-20, ifelse(is.na(par["K"]), 0, par["K"])-20)), Max=max(c(20, ifelse(is.na(par["K"]), 0, par["K"])+20)), Init=unname(ifelse(is.na(par["K"]), 0, par["K"])), row.names = "K"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["K1"]), 0, par["K1"])), Prior2=3, SDProp=0.5, Min=min(c(-20, ifelse(is.na(par["K1"]), 0, par["K1"])-20)), Max=max(c(20, ifelse(is.na(par["K1"]), 0, par["K1"])+20)), Init=unname(ifelse(is.na(par["K1"]), 0, par["K1"])), row.names = "K1"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["K2"]), 0, par["K2"])), Prior2=3, SDProp=0.5, Min=min(c(-20, ifelse(is.na(par["K2"]), 0, par["K2"])-20)), Max=max(c(20, ifelse(is.na(par["K2"]), 0, par["K2"])+20)), Init=unname(ifelse(is.na(par["K2"]), 0, par["K2"])), row.names = "K2"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=par["P_low"], Prior2=sdpx, SDProp=sdpx, Min=min(c(minx, par["P_low"]-5)), Max=max(c(maxx, par["P_low"]+5)), Init=par["P_low"], row.names ="P_low"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=par["P_high"], Prior2=sdpx, SDProp=sdpx, Min=min(c(minx, par["P_high"]-5)), Max=max(c(maxx, par["P_high"]+5)), Init=par["P_high"], row.names ="P_high"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["S_low"]), 0, par["S_low"])), Prior2=0.5, SDProp=0.5, Min=min(c(-2, ifelse(is.na(par["S_low"]), 0,par["S_low"])-2)), Max=max(c(2, ifelse(is.na(par["S_low"]), 0, par["S_low"])+2)), Init=unname(ifelse(is.na(par["S_low"]), 0, par["S_low"])), row.names = "S_low"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["S_high"]), 0, par["S_high"])), Prior2=0.5, SDProp=0.5, Min=min(c(-2, ifelse(is.na(par["S_high"]), 0,par["S_high"])-2)), Max=max(c(2, ifelse(is.na(par["S_high"]), 0, par["S_high"])+2)), Init=unname(ifelse(is.na(par["S_high"]), 0, par["S_high"])), row.names = "S_high"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["K1_low"]), 0, par["K1_low"])), Prior2=3, SDProp=0.5, Min=min(c(-20, ifelse(is.na(par["K1_low"]), 0, par["K1_low"])-20)), Max=max(c(20, ifelse(is.na(par["K1_low"]), 0, par["K1_low"])+20)), Init=unname(ifelse(is.na(par["K1_low"]), 0, par["K1_low"])), row.names = "K1_low"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["K2_low"]), 0, par["K2_low"])), Prior2=3, SDProp=0.5, Min=min(c(-20, ifelse(is.na(par["K2_low"]), 0, par["K2_low"])-20)), Max=max(c(20, ifelse(is.na(par["K2_low"]), 0, par["K2_low"])+20)), Init=unname(ifelse(is.na(par["K2_low"]), 0, par["K2_low"])), row.names = "K2_low"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["K1_high"]), 0, par["K1_high"])), Prior2=3, SDProp=0.5, Min=min(c(-20, ifelse(is.na(par["K1_high"]), 0, par["K1_high"])-20)), Max=max(c(20, ifelse(is.na(par["K1_high"]), 0, par["K1_high"])+20)), Init=unname(ifelse(is.na(par["K1_high"]), 0, par["K1_high"])), row.names = "K1_high"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["K2_high"]), 0, par["K2_high"])), Prior2=3, SDProp=0.5, Min=min(c(-20, ifelse(is.na(par["K2_high"]), 0, par["K2_high"])-20)), Max=max(c(20, ifelse(is.na(par["K2_high"]), 0, par["K2_high"])+20)), Init=unname(ifelse(is.na(par["K2_high"]), 0, par["K2_high"])), row.names = "K2_high"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["SL"]), 0, par["SL"])), Prior2=0.5, SDProp=0.5, Min=ifelse(equation == "flexit**", 0, min(c(-2, ifelse(is.na(par["SH"]), 0, par["SH"])-2))), Max=max(c(ifelse(equation == "flexit**", 10, 2), ifelse(is.na(par["SL"]), 0, par["SL"])+2)), Init=unname(ifelse(is.na(par["SL"]), 0, par["SL"])), row.names = "SL"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["SH"]), 0, par["SH"])), Prior2=0.5, SDProp=0.5, Min=ifelse(equation == "flexit**", 0, min(c(-2, ifelse(is.na(par["SH"]), 0, par["SH"])-2))), Max=max(c(ifelse(equation == "flexit**", 10, 2), ifelse(is.na(par["SH"]), 0, par["SH"])+2)), Init=unname(ifelse(is.na(par["SH"]), 0, par["SH"])), row.names = "SH"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["TransitionS"]), 0, par["TransitionS"])), Prior2=0.5, SDProp=0.5, Min=min(c(-5, ifelse(is.na(par["TransitionS"]), 0, par["TransitionS"])-5)), Max=max(c(5, ifelse(is.na(par["SH"]), 0, par["SH"])+5)), Init=unname(ifelse(is.na(par["TransitionS"]), 0, par["TransitionS"])), row.names = "TransitionS"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["SL_low"]), 0, par["SL_low"])), Prior2=0.5, SDProp=0.5, Min=min(c(ifelse(equation == "flexit**", 0, -2), ifelse(is.na(par["SL_low"]), 0, par["SL_low"])-2)), Max=max(c(ifelse(equation == "flexit**", 10, 2), ifelse(is.na(par["SL_low"]), 0, par["SL_low"])+2)), Init=unname(ifelse(is.na(par["SL_low"]), 0, par["SL_low"])), row.names = "SL_low"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["SH_low"]), 0, par["SH_low"])), Prior2=0.5, SDProp=0.5, Min=min(c(ifelse(equation == "flexit**", 0, -2), ifelse(is.na(par["SH_low"]), 0, par["SH_low"])-2)), Max=max(c(ifelse(equation == "flexit**", 10, 2), ifelse(is.na(par["SH_low"]), 0, par["SH_low"])+2)), Init=unname(ifelse(is.na(par["SH_low"]), 0, par["SH_low"])), row.names = "SH_low"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["TransitionS_low"]), 0, par["TransitionS_low"])), Prior2=0.5, SDProp=0.5, Min=min(c(-5, ifelse(is.na(par["TransitionS_low"]), 0, par["TransitionS_low"])-5)), Max=max(c(5, ifelse(is.na(par["TransitionS_low"]), 0, par["TransitionS_low"])+5)), Init=unname(ifelse(is.na(par["TransitionS_low"]), 0, par["TransitionS_low"])), row.names = "TransitionS_low"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["SL_high"]), 0, par["SL_high"])), Prior2=0.5, SDProp=0.5, Min=min(c(ifelse(equation == "flexit**", 0, -2), ifelse(is.na(par["SL_high"]), 0, par["SL_high"])-2)), Max=max(c(ifelse(equation == "flexit**", 10, 2), ifelse(is.na(par["SL_high"]), 0, par["SL_high"])+2)), Init=unname(ifelse(is.na(par["SL_high"]), 0, par["SL_high"])), row.names = "SL_high"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["SH_high"]), 0, par["SH_high"])), Prior2=0.5, SDProp=0.5, Min=min(c(ifelse(equation == "flexit**", 0, -2), ifelse(is.na(par["SH_high"]), 0, par["SH_high"])-2)), Max=max(c(ifelse(equation == "flexit**", 10, 2), ifelse(is.na(par["SH_high"]), 0, par["SH_high"])+2)), Init=unname(ifelse(is.na(par["SH_high"]), 0, par["SH_high"])), row.names = "SH_high"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=unname(ifelse(is.na(par["TransitionS_high"]), 0, par["TransitionS_high"])), Prior2=0.5, SDProp=0.5, Min=min(c(-5, ifelse(is.na(par["TransitionS_high"]), 0, par["TransitionS_high"])-5)), Max=max(c(5, ifelse(is.na(par["TransitionS_high"]), 0, par["TransitionS_high"])+5)), Init=unname(ifelse(is.na(par["TransitionS_high"]), 0, par["TransitionS_high"])), row.names = "TransitionS_high"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=ifelse(is.na(par["TRTL"]), 1, par["TRTL"]), Prior2=sdpx, SDProp=sdpx, Min=min(c(minx-10, par["TRTL"]-5)), Max=max(c(maxx+10, par["TRTL"]+5)), Init=par["TRTL"], row.names ="TRTL"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=ifelse(is.na(par["TRTH"]), 1, par["TRTH"]), Prior2=sdpx, SDProp=sdpx, Min=min(c(minx-10, par["TRTH"]-5)), Max=max(c(maxx+10, par["TRTH"]+5)), Init=par["TRTH"], row.names ="TRTH"))
  dpriors_norm <- rbind(dpriors_norm, 
                        data.frame(Density="dnorm", Prior1=ifelse(is.na(par["SRTRT"]), 0, par["SRTRT"]), Prior2=2, SDProp=1, Min=min(c(-10, par["SRTRT"]-5)), Max=max(c(10, par["SRTRT"]+5)), Init=par["SRTRT"], row.names ="SRTRT"))
  
  dpriors_unif <- dpriors_norm
  dpriors_unif[, "Density"] <- "dunif"
  dpriors_unif[, "Prior1"] <- dpriors_norm[, "Min"]
  dpriors_unif[, "Prior2"] <- dpriors_norm[, "Max"]
  
  parameters <- NULL
  for (np in seq_along(par)) {
    if (default[np] == "dnorm") {
      parameters <- rbind(parameters, dpriors_norm[names(par[np]), ])
    } else {
      parameters <- rbind(parameters, dpriors_unif[names(par[np]), ])
    }
  }
  
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
