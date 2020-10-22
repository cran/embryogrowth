#' HatchingSuccess.MHmcmc_p generates set of parameters to be used with HatchingSuccess.MHmcmc()
#' @title Generates set of parameters to be used with HatchingSuccess.MHmcmc()
#' @author Marc Girondot
#' @return A matrix with the parameters
#' @param result An object obtained after a HatchingSuccess.fit() fit
#' @param parameters A set of parameters. Replace the one from result
#' @param fixed.parameters A set of fixed parameters. Replace the one from result
#' @param accept If TRUE, the script does not wait user information
#' @description Interactive script used to generate set of parameters to be used with HatchingSuccess.MHmcmc().\cr
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
#' g <- HatchingSuccess.fit(par=par, data=totalIncubation_Cc)
#' pMCMC <- HatchingSuccess.MHmcmc_p(g, accept=TRUE)
#' mcmc <- HatchingSuccesss.MHmcmc(result=g, parameters = pMCMC, 
#'                   adaptive=TRUE, n.iter=100000, trace=1000)
#' }
#' @export

# Algo Metropolis-Hastings
# ------------------------

HatchingSuccess.MHmcmc_p<-function(result=NULL, parameters=NULL, fixed.parameters=NULL, 
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

# 7/2/2014, ajout de la nouvelle version des parametres
# #' par <- c(S.low=0.5, S.high=0.3,  P.low=25, deltaP=10, MaxHS=logit(0.8)) 
  
  P.low <- abs(allpar["P.low"])
  deltaP <- abs(allpar["deltaP"])
  P.high <- abs(allpar["P.high"])
  
  S.low <- allpar["S.low"]
  S.high <- allpar["S.high"]
  
  K1.low <- allpar["K1.low"]
  K1.high <- allpar["K1.high"]
  K2.low <- allpar["K2.low"]
  K2.high <- allpar["K2.high"]  
  
  if (is.na(P.low)) P.low <- P.high - deltaP
  if (is.na(P.high)) P.high <- P.low + deltaP
  if (is.na(deltaP)) deltaP <- P.high - P.low
  
  
  # S.low
  S.low <- c("dunif", 0, max(allpar["S.low"]*2, 10), 2, 0, max(allpar["S.low"]*2, 10), par["S.low"])
  
  # S.high
  S.high <- c("dunif", 0, max(allpar["S.high"]*2, 10), 2, 0, max(par["S.high"]*2, 10), par["S.high"])
  
  # P.low
  P.low <- c("dunif", 0, max(P.low * 2, 30), 2, 0, max(P.low * 2, 30), P.low)

  # P.high
  P.high <- c("dunif", 0, max(P.high * 2, 40), 2, 0, max(P.high * 2, 40), P.high)
  
  # deltaP
  deltaP <- c("dunif", 0, max(deltaP * 2, 10), 2, 0, max(deltaP * 2, 10), deltaP)
  
  # MaxHS
  MaxHS <- c("dunif", 0, 1, 2, 0, 1, par["MaxHS"])
  
  K1.low <- c("dunif", min(K1.low * 2, -10), max(K1.low * 2, 10), 2, min(K1.low * 2, -10), max(K1.low * 2, 10), K1.low)
  K1.high <- c("dunif", min(K1.high * 2, -10), max(K1.high * 2, 10), 2, min(K1.high * 2, -10), max(K1.high * 2, 10), K1.high)
  K2.low <- c("dunif", min(K2.low * 2, -10), max(K2.low * 2, 10), 2, min(K2.low * 2, -10), max(K2.low * 2, 10), K2.low)
  K2.high <- c("dunif", min(K2.high * 2, -10), max(K2.high * 2, 10), 2, min(K2.high * 2, -10), max(K2.high * 2, 10), K2.high)
  
  
priors <- list(S.low=S.low, S.high=S.high, P.low=P.low, P.high=P.high, deltaP=deltaP, MaxHS=MaxHS, 
               K1.low=K1.low, K1.high=K1.high, K2.low=K2.low, K2.high=K2.high)

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
