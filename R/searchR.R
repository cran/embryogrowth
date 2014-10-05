#' searchR fits the parameters that best represent nest incubation data.
#' @title Fit the parameters that best represent nest incubation data.
#' @author Marc Girondot
#' @return A result object
#' @param parameters A set of parameters used as initial point for searching
#' @param fixed.parameters A set of parameters that will not be changed
#' @param temperatures Timeseries of temperatures after formated using FormatNests()
#' @param derivate Function used to fit embryo growth: dydt.Gompertz, dydt.exponential or dydt.linear
#' @param test Mean and SD of size of hatchlings
#' @param M0 Measure of hatchling size or mass proxi at laying date
#' @param method Method uses for searching. Can be any method from optim function
#' @param maxiter After maxiter iteration, the value of parameters is displayed but it continues if convergence is not acheived
#' @param saveAtMaxiter If True, each time number of interation reach maxiter, current data are saved in file with filename name
#' @param fileName The intermediate results are saved in file with fileName.Rdata name
#' @param weight A named vector of the weight for each nest for likelihood estimation
#' @param hessian If TRUE, the hessian matrix is estimated and the SE of parameters estimated.
#' @param parallel If true, try to use several cores using parallel computing. Must be FALSE in Windows.
#' @description Fit the parameters that best represent data.
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' # K for Gompertz must be set as fixed parameter or being a constant K  
#' # or relative to the hatchling size rK
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' resultNest_4p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92), method = "BFGS", maxiter = 200)
#' data(resultNest_4p)
#' plot(resultNest_4p, xlim=c(0,70), ylimT=c(22, 32), ylimS=c(0,45), series=1)
#' x <- structure(c(115.758929130522, 428.649022170996, 503.687251738993, 
#' 12.2621455821612, 306.308841227278, 116.35048615105), .Names = c("DHA", 
#' "DHH", "DHL", "DT", "T12L", "Rho25"))
#' pfixed <- c(rK=2.093313)
#' resultNest_6p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92), method = "BFGS", maxiter = 200)
#' data(resultNest_6p)
#' pMCMC <- embryogrowth_MHmcmc_p(resultNest_6p, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' result_mcmc_6p <- embryogrowth_MHmcmc(result=resultNest_6p,  
#' 	parametersMCMC=pMCMC, n.iter=10000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(result_mcmc_6p)
#' # compare_AIC() is a function from the package "phenology"
#' outputAIC<-compare_AIC(list(test1=resultNest_4p, test2=resultNest_6p))
#' ############ with new parametrization
#' data(resultNest_4p)
#' x0 <- resultNest_4p$par
#' t <- hist(resultNest_4p, plot=FALSE)
#' temperatures <- seq(from=floor(range(t$temperatures)[1]+273.15-1), 
#' to=floor(range(t$temperatures)[2]+273.15+1), length.out=7)
#' newx <- embryogrowth:::.SSM(temperatures, x0)[[1]]*1E5
#' names(newx) <- temperatures
#' pfixed <- c(rK=2.093313)
#' resultNest_newp <- searchR(parameters=newx, fixed.parameters=pfixed,
#'  temperatures=formated, derivate=dydt.Gompertz, M0=1.7,
#'  test=c(Mean=39.33, SD=1.92), method = "BFGS", maxiter = 200)
#' plotR_hist(resultNest_newp, ylim=c(0,0.3), xlimR=c(23, 34), ylimH=c(0, 0.3))
#' outputAIC<-compare_AIC(list(test4p=resultNest_4p, test6p=resultNest_6p, 
#'  testAnchor=resultNest_newp))
#' }
#' @export


searchR <-
function(parameters=stop('Initial set of parameters must be provided'), 
	fixed.parameters=NULL, temperatures=stop('Temperature data must be provided !'), 
	derivate=dydt.Gompertz, test=c(Mean=39.33, SD=1.92), 
	M0=1.7, method="BFGS", maxiter=200, saveAtMaxiter=FALSE, fileName="intermediate", 
  weight=NULL, hessian=TRUE, parallel=(.Platform$OS.type=="unix")) {
	
	if (parallel & (.Platform$OS.type!="unix")) {
		parallel <- FALSE
		warning("Parallel computing is not available for Windows")
	}
  
	if (!requireNamespace("numDeriv", quietly = TRUE)) {
	  warning("numDeriv package is absent; less accurate fitting method will be used and SE of parameters cannot be calculated")
	  numDeriv <- FALSE
    hessian <- FALSE
	} else {
	  numDeriv <- TRUE
	}
	
	
NbTS <- temperatures$IndiceT[3]
	
if (is.null(weight)) {
	par <- rep(1, NbTS)
	names(par) <- names(temperatures)[1:NbTS]
} else {

	if (any(is.na(weight))) {
		par <- rep(1, NbTS)
		names(par) <- names(temperatures)[1:NbTS]
	} else {

	if (is.list(weight)) weight <- weight$weight

	if (length(setdiff(names(weight), names(temperatures)[1:NbTS]))==0) {
		par <- weight
	} else {
		print("Check the weights")
		return(invisible())
	}
	}
}

weight <- par

# test si tous sont là
if (length(setdiff(names(temperatures)[1:temperatures$IndiceT[3]], names(weight)))!=0) {
	print("The weight parameter must define weight for each nest.")
	print(paste("check", setdiff(names(temperatures)[1:temperatures$IndiceT[3]], names(weight)), "nests"))
	return(invisible())	
}	


##########################################################
# Données de base de Gompertz
##########################################################

if (is.numeric(test)) {
	testuse<-data.frame(Mean=rep(test["Mean"], temperatures[["IndiceT"]][3]), SD=rep(test["SD"], temperatures[["IndiceT"]][3]), row.names=names(temperatures[1:temperatures$IndiceT["NbTS"]]))
} else {
	testuse<-test
}

for (j in 1:temperatures[["IndiceT"]][3]) temperatures[[names(temperatures)[j]]][1, "Mass"] <- M0



# Un paramètre ne peut pas être indiqué en fixe et en fité - 22/7/2012	
# test faux, corrigé le 19/2/2013
	if (length(intersect(names(parameters), names(fixed.parameters)))!=0) {
		print("A parameter cannot be fixed and fitted at the same time !")
		return(invisible())
	}

ptec <- list(temperatures=temperatures, derivate=derivate, weight=weight,
		testuse=testuse, M0=M0, fixed.parameters=fixed.parameters, parallel=parallel)

repeat {

# 30/8/2014 - Je retire optimx
  if (numDeriv) {
	result=optim(parameters, .fonctionfit, pt=ptec,
               gr=.gradientRichardson, method=method, 
               control=list(trace=2, REPORT=1, maxit=maxiter), hessian=FALSE)
} else {
  result=optim(parameters, .fonctionfit, pt=ptec,
               method=method, control=list(trace=2, REPORT=1, maxit=maxiter), hessian=FALSE)
}
	if (result$convergence==0) break
	parameters<-result$par
	print("Convergence is not achieved. Optimization continues !")
	print(parameters)
  if (saveAtMaxiter) save(parameters, file=paste0(fileName, ".RData"))
}

print(result$par)

if (hessian) {

	mathessian <- try(hessian(.fonctionfit, result$par, method="Richardson", pt=ptec), silent=TRUE)

	if (inherits(mathessian, "try-error")) {
			res<-rep(NA, length(parameters))
			result$hessian <- NA
			print("SE of parameters cannot be estimated.")
			print("Probably the model is badly fitted. Try using searchR() before.")

	} else {
	

	result$hessian <- mathessian

	inversemathessian <- try(solve(mathessian), silent=TRUE)
	
	if (inherits(inversemathessian, "try-error")) {
		res <- -1
	} else {
		res <- diag(inversemathessian)
	}

	# Je gère plus correctement les erreurs - 17/7/2012

	neg=any(res<0)
	if (!neg) {
		res=sqrt(res)
	} else {
		res<-rep(NA, length(parameters))
		print("SE of parameters cannot be estimated.")
		if (inherits(inversemathessian, "try-error")) {
			print("Probably the model is badly fitted. Try other initial points.")
		} else {
			print("Probably flat likelihood is observed around some parameters.")
			print("Try using embryogrowth_MHmcmc() function to get the SE of parameters.")
		}
	}
}
	
} else {
# pas de hessian donc pas de SE
	res<-rep(NA, length(parameters))
}


names(res)=names(parameters)

result$SE <- res 

result$data <- temperatures

# Avant *5. Correction du 17/7/2012
result$AIC <- 2*result$value+2*length(parameters)

result$test <- testuse
result$derivate <- derivate
result$M0 <- M0
# Je stocke aussi les paramètres fixé-16/7/2012
result$fixed.parameters <- fixed.parameters
# 29/1/2014
result$weight <- weight

growlnotify("End of optimization !")

class(result) <- "NestsResult"

return(result)

}
