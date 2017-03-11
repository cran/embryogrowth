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
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' resultNest_4p_SSM4p <- searchR(parameters=x, fixed.parameters=pfixed,  
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7,  
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p_SSM4p)
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p_SSM4p, accept=TRUE)
#' }
#' @export

# Algo Metropolis-Hastings
# ------------------------

TRN_MHmcmc_p<-function(result=NULL, parameters=NULL, fixed.parameters=NULL, 
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
if (all(names(allpar)!="Rho25") & all(names(allpar)!="Peak") & all(names(allpar)!="k")) {
	priors <- list()
	for(i in 1:length(par)) {
		pr <- c("dunif", 0.001, 20, 1, 0.001, 20, par[i])
		priors <- c(priors, list(pr))
	}
	names(priors) <- names(par)

} else {
  
# k
  if (!is.na(par["k"])) par["k"] <- abs(par["k"])
  k <- c("dunif", 0, ifelse(is.na(par["k"]), 500, par["k"]*2), 2, 0, ifelse(is.na(par["k"]), 500, par["k"]*2), ifelse(is.na(par["k"]), 400, par["k"]))
  
# lambda
  if (!is.na(par["lambda"])) par["lambda"] <- abs(par["lambda"])
  lambda <- c("dunif", 0, ifelse(is.na(par["lambda"]), 500, par["lambda"]*2), 2, 0, ifelse(is.na(par["lambda"]), 500, par["lambda"]*2), ifelse(is.na(par["lambda"]), 400, par["lambda"]))
  
# scale
  if (!is.na(par["scale"])) par["scale"] <- abs(par["scale"])
  scale <- c("dunif", 0, ifelse(is.na(par["scale"]), 500, par["scale"]*2), 2, 0, ifelse(is.na(par["scale"]), 500, par["scale"]*2), ifelse(is.na(par["scale"]), 400, par["scale"]))

  # k_L
  if (!is.na(par["k_L"])) par["k_L"] <- abs(par["k_L"])
  k_L <- c("dunif", 0, ifelse(is.na(par["k_L"]), 500, par["k_L"]*2), 2, 0, ifelse(is.na(par["k_L"]), 500, par["k_L"]*2), ifelse(is.na(par["k_L"]), 400, par["k_L"]))
  
  # lambda_L
  if (!is.na(par["lambda_L"])) par["lambda_L"] <- abs(par["lambda_L"])
  lambda_L <- c("dunif", 0, ifelse(is.na(par["lambda_L"]), 500, par["lambda_L"]*2), 2, 0, ifelse(is.na(par["lambda_L"]), 500, par["lambda_L"]*2), ifelse(is.na(par["lambda_L"]), 400, par["lambda_L"]))
  
  # scale_L
  if (!is.na(par["scale_L"])) par["scale_L"] <- abs(par["scale_L"])
  scale_L <- c("dunif", 0, ifelse(is.na(par["scale_L"]), 500, par["scale_L"]*2), 2, 0, ifelse(is.na(par["scale_L"]), 500, par["scale_L"]*2), ifelse(is.na(par["scale_L"]), 400, par["scale_L"]))
  
# Peak
  Peak <- c("dnorm", 32, 2, 1, 20, 40, ifelse(is.na(par["Peak"]), 32, par["Peak"]))
  
# Length
  Length <- c("dnorm", 5, 2, 1, 0, 20, ifelse(is.na(par["Length"]), 10, par["Length"]))

# LengthB
  LengthB <- c("dnorm", 5, 2, 1, 0, 20, ifelse(is.na(par["LengthB"]), 10, par["LengthB"]))

# LengthE
  LengthE <- c("dnorm", 5, 2, 1, 0, 20, ifelse(is.na(par["LengthE"]), 10, par["LengthE"]))
  
# Max
  Max <- c("dnorm", 5, 2, 1, 0, 20, ifelse(is.na(par["Max"]), 10, par["Max"]))
  
# sd
  sd <- c("dnorm", 5, 2, 1, 0, 20, ifelse(is.na(par["sd"]), 10, par["sd"]))
  
# sdH
  sdH <- c("dnorm", 5, 2, 1, 0, 20, ifelse(is.na(par["sdH"]), 10, par["sdH"]))
  
# sdL
  sdL <- c("dnorm", 5, 2, 1, 0, 20, ifelse(is.na(par["sdL"]), 10, par["sdL"]))

# Scale
  Scale <- c("dnorm", 3, 5, 1, 0, 100, ifelse(is.na(par["Scale"]), 3, par["Scale"]))
  
# "epsilon"
  
  epsilon <- c("dunif", -10, 10, 0.5, -10, 10, 0)

# "DHA"

DHA <- c("dunif", 1, 500, 2, 1, 500, ifelse(is.na(par["DHA"]), 200, par["DHA"]))
if (as.numeric(DHA[7])>400) {
	DHA <- c("dunif", 1, as.numeric(DHA[7])+200, 2, 1, as.numeric(DHA[7])+200, as.numeric(DHA[7]))
}

# "DHL"

DHL <- c("dunif", 1, 500, 2, 1, 500, ifelse(is.na(par["DHL"]), 370, par["DHL"]))
if (as.numeric(DHL[7])>400) {
	DHL <- c("dunif", 1, as.numeric(DHL[7])+200, 2, 1, as.numeric(DHL[7])+200, as.numeric(DHL[7]))
}


# "DHH"

DHH <- c("dunif", 1, 500, 2, 1, 500, ifelse(is.na(par["DHH"]), 400, par["DHH"]))
if (as.numeric(DHH[7])>400) {
	DHH <- c("dunif", 1, as.numeric(DHH[7])+200, 2, 1, as.numeric(DHH[7])+200, as.numeric(DHH[7]))
}


# "T12L"

T12L <- c("dnorm", 30+273.15, 7, 1, 20+273.15, 40+273.15, ifelse(is.na(par["T12L"]), 30+273.15, par["T12L"]))

# "T12H"

T12H <- c("dnorm", 30+273.15, 7, 1, 20+273.15, 40+273.15, ifelse(is.na(par["T12H"]), 30+273.15, par["T12H"]))

# "DT"

DT <- c("dnorm", 15, 10, 1, 0, 100, ifelse(is.na(par["DT"]), 10, par["DT"]))

# "Rho25"

Rho25 <- c("dunif", 1, 500, 2, 1, 500, ifelse(is.na(par["Rho25"]), 100, par["Rho25"]))
if (as.numeric(Rho25[7])>400) {
	Rho25 <- c("dunif", 1, as.numeric(Rho25[7])+200, 2, 1, as.numeric(Rho25[7])+200, as.numeric(Rho25[7]))
}

# "DHA_L"

DHA_L <- c("dunif", 1, 500, 2, 1, 500, ifelse(is.na(par["DHA_L"]), 200, par["DHA_L"]))
if (as.numeric(DHA_L[7])>400) {
	DHA_L <- c("dunif", 1, as.numeric(DHA[7])+200, 2, 1, as.numeric(DHA_L[7])+200, as.numeric(DHA_L[7]))
}

# "DHL_L"

DHL_L <- c("dunif", 1, 500, 2, 1, 500, ifelse(is.na(par["DHL_L"]), 370, par["DHL_L"]))
if (as.numeric(DHL_L[7])>400) {
	DHL_L <- c("dunif", 1, as.numeric(DHL_L[7])+200, 2, 1, as.numeric(DHL_L[7])+200, as.numeric(DHL_L[7]))
}


# "DHH_L"

DHH_L <- c("dunif", 1, 500, 2, 1, 500, ifelse(is.na(par["DHH_L"]), 400, par["DHH_L"]))
if (as.numeric(DHH_L[7])>400) {
	DHH_L <- c("dunif", 1, as.numeric(DHH_L[7])+200, 2, 1, as.numeric(DHH_L[7])+200, as.numeric(DHH_L[7]))
}


# "T12L_L"

T12L_L <- c("dnorm", 30+273.15, 7, 1, 20+273.15, 40+273.15, ifelse(is.na(par["T12L_L"]), 30+273.15, par["T12L_L"]))

# "T12H_L"

T12H_L <- c("dnorm", 30+273.15, 7, 1, 20+273.15, 40+273.15, ifelse(is.na(par["T12H_L"]), 30+273.15, par["T12H_L"]))

# "DT"

DT_L <- c("dnorm", 15, 10, 1, 0, 100, ifelse(is.na(par["DT_L"]), 10, par["DT_L"]))

# "Rho25"

Rho25_L <- c("dunif", 1, 500, 2, 1, 500, ifelse(is.na(par["Rho25_L"]), 100, par["Rho25_L"]))
if (as.numeric(Rho25_L[7])>400) {
	Rho25_L <- c("dunif", 1, as.numeric(Rho25_L[7])+200, 2, 1, as.numeric(Rho25_L[7])+200, as.numeric(Rho25_L[7]))
}



# "transition_P"

transition_P <- c("dunif", 0, 40, 2, 0, 40, ifelse(is.na(par["transition_P"]), 20, par["transition_P"]))

# "transition_S"

transition_S <- c("dunif", -100, 100, 2, -100, 100, ifelse(is.na(par["transition_S"]), 4, par["transition_S"]))

priors <- list(DHA=DHA, DHL=DHL, DHH=DHH, T12L=T12L, T12H=T12H, DT=DT, Rho25=Rho25, 
DHA_L=DHA_L, DHL_L=DHL_L, DHH_L=DHH_L, T12L_L=T12L_L, T12H_L=T12H_L, DT_L=DT_L, Rho25_L=Rho25_L,
transition_P=transition_P, transition_S=transition_S, epsilon=epsilon, 
Peak=Peak, Scale=Scale, sdL=sdL, sdH=sdH, sd=sd, Length=Length, LengthB=LengthB, 
LengthE=LengthE, Max=Max, 
k=k, scale=scale, lambda=lambda, 
k_L=k_L, scale_L=scale_L, lambda_L=lambda_L)

}

prencours <- NULL

for (i in 1:length(par)) {

	prencours <- c(prencours, priors[[names(par)[i]]])
}

parametersMCMC <- matrix(prencours, ncol=7, byrow=T)
colnames(parametersMCMC) <- c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", "Init")
rownames(parametersMCMC)<-names(par)

parametersMCMC <- as.data.frame(parametersMCMC, stringsAsFactors = FALSE)

for (i in 2:7)
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
