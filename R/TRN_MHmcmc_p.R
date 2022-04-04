#' TRN_MHmcmc_p generates set of parameters to be used with GRTRN_MHmcmc() or STRN_MHmcmc()
#' @title Generates set of parameters to be used with GRTRN_MHmcmc() or STRN_MHmcmc()
#' @author Marc Girondot
#' @return A matrix with the parameters
#' @param result An object obtained after a SearchR fit
#' @param parameters A set of parameters. Replace the one from result
#' @param fixed.parameters A set of fixed parameters. Replace the one from result
#' @param accept If TRUE, the script does not wait user information
#' @description Interactive script used to generate set of parameters to be used with GRTRN_MHmcmc() or STRN_MHmcmc().\cr
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
#' x <- structure(c(105.966881676793, 613.944134764125, 306.449533440186, 
#'                 118.193882815108), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' resultNest_4p_SSM <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, integral=integral.Gompertz, M0=M0, 
#' 	hatchling.metric=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p_SSM)
#' plot(resultNest_4p_SSM, xlim=c(0,70), ylimT=c(22, 32), ylimS=c(0,45), series=1, 
#' embryo.stages="Caretta caretta.SCL")
#' ############################################################################
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p_SSM, accept=TRUE)
#' }
#' @export

# Algo Metropolis-Hastings
# ------------------------

TRN_MHmcmc_p <- function(result=NULL, parameters=NULL, fixed.parameters=NULL, 
                         accept=FALSE) {
  
  # result=NULL; parameters=NULL; fixed.parameters=NULL; accept=TRUE
  
  # d'abord je sors les parametres a utiliser
  
  if (is.null(result) & is.null(parameters)) {
    stop("Or result or parameters must be provided")
  }
  
  # 26/4/2015
  if (is.null(parameters)) parameters <- result$par
  if (is.null(fixed.parameters)) fixed.parameters <- result$fixed.parameters
  
  par <- parameters
  allpar <- c(parameters, fixed.parameters)
  
  # 2/4/2022, Je simplifie tout
  priors_limits <- data.frame(DHA=c(-100, 500), 
                              DHL=c(-100, 500), 
                              DHH=c(-100, 500), 
                              T12L=c(-100, 500), 
                              T12H=c(-100, 500), 
                              DT=c(0, 500), 
                              Rho25=c(0, 500), 
                              DHA_L=c(-100, 500), 
                              DHL_L=c(-100, 500), 
                              DHH_L=c(-100, 500), 
                              T12L_L=c(-100, 500), 
                              T12H_L=c(-100, 500), 
                              DT_L=c(0, 500), 
                              Rho25_L=c(0, 500),
                              transition_P=c(-100, 500), 
                              transition_S=c(-100, 500), 
                              epsilon=c(0, 500), 
                              SD=c(0, 500), 
                              Peak=c(-100, 500), 
                              Scale=c(0, 500), 
                              sdL=c(0, 500), 
                              sdH=c(0, 500), 
                              sd=c(0, 500), 
                              Length=c(0, 500), 
                              LengthB=c(0, 500), 
                              LengthE=c(0, 500), 
                              Max=c(0, 500), 
                              k=c(0, 500), 
                              scale=c(0, 500), 
                              lambda=c(0, 500), 
                              k_L=c(0, 500), 
                              scale_L=c(0, 500), 
                              lambda_L=c(0, 500))
  parameters <- NULL
    for(i in 1:length(par)) {
      minp <- priors_limits[1, names(par[i])]
      maxp <- priors_limits[2, names(par[i])]
      parameters <- rbind(parameters, 
                      data.frame(Density="dunif", 
                                 Prior1=min(c(par[i]-abs(par[i])*0.5, minp)), 
                                 Prior2=max(c(par[i]+abs(par[i])*0.5, maxp)), 
                                 SDProp=1, 
                                 Min=min(c(par[i]-abs(par[i])*0.5, minp)), 
                                 Max=max(c(par[i]+abs(par[i])*0.5, maxp)), 
                                 Init=par[i], 
                                 row.names = names(par[i])))
    }
    
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
            
            cat(paste("Distribution of the prior, Minimum (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
            f<-scan(nmax=1, quiet=TRUE, what=character())
            if (length(f)!=0) parameters[variable, "Prior1"] <- f
            cat(paste("Distribution of the prior, Maximum (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
            f<-scan(nmax=1, quiet=TRUE, what=character())
            if (length(f)!=0) parameters[variable, "Prior2"] <- f
            
          } else {
            
            if (density == "dnorm") {
              
              cat(paste("Distribution of the prior, Mean (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
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
