#' searchR fits the parameters that best represent nest incubation data.
#' @title Fit the parameters that best represent nest incubation data.
#' @author Marc Girondot
#' @return A result object
#' @param parameters A set of parameters used as initial point for searching
#' @param fixed.parameters A set of parameters that will not be changed
#' @param temperatures Timeseries of temperatures after formated using FormatNests()
#' @param derivate Function used to fit embryo growth: dydt.Gompertz, dydt.exponential or dydt.linear
#' @param test A vector with Mean and SD of size of hatchlings, ex. test=c(Mean=39, SD=3)
#' @param M0 Measure of hatchling size or mass proxi at laying date
#' @param method Method uses for searching. Can be any method from optim function
#' @param maxiter After maxiter iteration, the value of parameters is displayed but it continues if convergence is not acheived
#' @param saveAtMaxiter If True, each time number of interation reach maxiter, current data are saved in file with filename name
#' @param fileName The intermediate results are saved in file with fileName.Rdata name
#' @param weight A named vector of the weight for each nest for likelihood estimation
#' @param hessian If TRUE, the hessian matrix is estimated and the SE of parameters estimated.
#' @description Fit the parameters that best represent data.\cr
#' test can be a list with two elements Mean and SD and each element is a nammed vector with the nest name.
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
#' x <- structure(c(118.431040984352, 498.205702157603, 306.056280989839, 
#' 118.189669472381), .Names = c("DHA", "DHH", "T12H", "Rho25"))
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
#' pMCMC <- TRN_MHmcmc_p(resultNest_6p, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' result_mcmc_6p <- GRTRN_MHmcmc(result=resultNest_6p,  
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
	fixed.parameters=NULL, temperatures=stop('Formated temperature must be provided !'), 
	derivate=dydt.Gompertz, test=c(Mean=39.33, SD=1.92), 
	M0=1.7, method="BFGS", maxiter=200, saveAtMaxiter=FALSE, fileName="intermediate", 
  weight=NULL, hessian=TRUE) {
  
  # parameters <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25")); fixed.parameters <- c(rK=2.093313)
  # temperatures <- formated;derivate <- dydt.Gompertz;M0 <- 1.7
  # test=c(Mean=39.33, SD=1.92); method = "BFGS"; maxiter = 200;saveAtMaxiter=FALSE;fileName="intermediate"
	# weight=NULL; hessian=TRUE
  
	if (!requireNamespace("numDeriv", quietly = TRUE)) {
	  warning("numDeriv package is absent; less accurate fitting method will be used and SE of parameters cannot be calculated")
	  numDeriv <- FALSE
    hessian <- FALSE
	} else {
	  numDeriv <- TRUE
	}
	
  # dans temperatures il faut que je rajoute une colonne avec les indices de températures en K
  
	
NbTS <- temperatures$IndiceT[3]

# si j'ai weight dans les data formatées et pas en paramètres, je les prends
if (is.null(weight) & !is.null(temperatures$weight)) {
  weight <- temperatures$weight
}

# si j'ai pas de weight, je mets tout à 1
if (is.null(weight)) {
  weight <- rep(1, NbTS)
	names(weight) <- names(temperatures)[1:NbTS]
}

# si c'est une liste, je prends l'élément weight
if (is.list(weight)) weight <- weight$weight

if (!setequal(names(weight), names(temperatures)[1:NbTS])) {
  print("The weight parameter must define weight for each nest.")
  print(paste("check", setdiff(names(temperatures)[1:temperatures$IndiceT[3]], names(weight)), "nests"))
  return(invisible())	
}

##########################################################
# Données de base de Gompertz
##########################################################

if (is.numeric(test)) {
	testuse<-data.frame(Mean=rep(test["Mean"], NbTS), SD=rep(test["SD"], NbTS), row.names=names(temperatures[1:NbTS]))
} else {
	testuse<-test
}

# 25/2/2015
for (j in 1:NbTS) temperatures[[j]][1, "Mass"] <- M0

# Un paramètre ne peut pas être indiqué en fixe et en fité - 22/7/2012	
# test faux, corrigé le 19/2/2013
	if (length(intersect(names(parameters), names(fixed.parameters)))!=0) {
		print("A parameter cannot be fixed and fitted at the same time !")
		return(invisible())
	}

repeat {

# 30/8/2014 - Je retire optimx
  if (numDeriv) {
	result=optim(parameters, .fonctionfit, temperatures=temperatures, 
	             derivate=derivate, weight=weight,
	             test=testuse, M0=M0, fixed.parameters=fixed.parameters,
               gr=.gradientRichardson, method=method, 
               control=list(trace=2, REPORT=1, maxit=maxiter), hessian=FALSE)
} else {
  result=optim(parameters, .fonctionfit, temperatures=temperatures, 
               derivate=derivate, weight=weight,
               test=testuse, M0=M0, fixed.parameters=fixed.parameters,
               method=method, control=list(trace=2, REPORT=1, maxit=maxiter), hessian=FALSE)
}
	if (result$convergence==0) break
	parameters<-result$par
	print("Convergence is not achieved. Optimization continues !")
	print(dput(parameters))
 if (saveAtMaxiter) save(parameters, file=paste0(fileName, ".RData"))
}

print(result$par)

if (hessian) {

	mathessian <- try(hessian(.fonctionfit, result$par, method="Richardson", temperatures=temperatures, 
	                          derivate=derivate, weight=weight,
	                          test=testuse, M0=M0, fixed.parameters=fixed.parameters), silent=TRUE)

	if (inherits(mathessian, "try-error")) {
			res<-rep(NA, length(parameters))
			result$hessian <- NA
			print("SE of parameters cannot be estimated.")
			print("Probably the model is badly fitted. Try using searchR() again or use GRTRN_MHmcmc().")

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
			print("Try using GRTRN_MHmcmc() function to get the SE of parameters.")
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
