#' STRN estimates the parameters that best describe the sexualisation thermal reaction norm within the TSP
#' @title Estimate the parameters that best describe the sexualisation thermal reaction norm within the TSP
#' @author Marc Girondot
#' @return The list with object return by optim() 
#' @param Initial_STRN Values for initial model of Sexualisation Thermal Reaction Norm
#' @param fixed.parameters Value for Sexualisation Thermal Reaction Norm model that will not be changed
#' @param EmbryoGrowthTRN The Embryo Growth Thermal Reaction Norm obtained with searchR()
#' @param tsd The model used to predict sex ratio, obtained from tsd()
#' @param Sexed The number of sexed embryos with names identifying timeseries
#' @param Males The number of males embryos with names identifying timeseries
#' @param Females The number of females embryos with names identifying timeseries
#' @param Temperatures The temperature from out of info.nests to be used
#' @param SE Should standard error of parameters and Hessian matrix be estimated ? TRUE or FALSE
#' @param method Methods to be used with optimx
#' @param itnmax Maximum number of iterations for each method
#' @param parallel Should parallel computing for info.nests() be used
#' @param control List for control parameters for optimx
#' @description Estimate the parameters that best describe the sexualisation thermal reaction norm within the TSP.\cr
#' The Temperatures parameter is a character string which can be:\cr
#' \itemize{
#'   \item \code{TimeWeighted.temperature.mean}
#'   \item \code{TSP.TimeWeighted.temperature.mean}
#'   \item \code{TSP.MassWeighted.temperature.mean}
#'   \item \code{TSP.STRNWeighted.temperature.mean}
#'   \item \code{TSP.MassWeighted.STRNWeighted.temperature.mean}
#'   \item \code{MiddleThird.TimeWeighted.temperature.mean}
#'   }
#' If information for sex is not known for some timeseries, set NA for Sexed.\cr
#' Sexed, Males and Females must be vectors with names. The names must be the same as 
#' the names of timeseries of temperatures in EmbryoGrowthTRN.\cr
#' Only two of these 3 parameters are required: Males, Females and Sexed\cr
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' MedIncubation_Cc <- subset(DatabaseTSD, Species=="Caretta caretta" & 
#' RMU=="Mediterranean" & Sexed!=0)
#' Med_Cc <- with(MedIncubation_Cc, tsd(males=Males, females=Females, 
#'  temperatures=Incubation.temperature, par=c(P=29.5, S=-0.01)))
#' plot(Med_Cc, xlim=c(25, 35))
#' males <- c(7, 0, 0, 0, 0, 5, 6, 3, 5, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' names(males) <- rev(rev(names(resultNest_4p_SSM4p$data))[-(1:2)])
#' sexed <- rep(10, length(males))
#' names(sexed) <- rev(rev(names(resultNest_4p_SSM4p$data))[-(1:2)])
#' Initial_STRN <- resultNest_4p_SSM4p$par[c("DHA", "DHH", "T12H")]
#' Initial_STRN <- structure(c(582.567096666926, 2194.0806711639, 3475.28414940385), 
#'                           .Names = c("DHA", "DHH", "T12H"))
#' fp <- c(Rho25=100)
#' fitSTRN <- STRN(Initial_STRN=Initial_STRN, 
#'                 EmbryoGrowthTRN=resultNest_4p_SSM4p, tsd=Med_Cc, 
#'                 Sexed=sexed, Males=males, 
#'                 fixed.parameters=fp, 
#'                 SE=TRUE, 
#'                 Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean")
#' plotR(fitSTRN, curves ="ML quantiles", ylim=c(0,2))
#' CTE <- info.nests(NestsResult=resultNest_4p_SSM4p, 
#'                   SexualisationTRN=fitSTRN,
#'                   SexualisationTRN.CI="Hessian",
#'                   CI="Hessian", 
#'                   replicate.CI=100, 
#'                   progress=TRUE, 
#'                   warnings=TRUE, 
#'                   out="summary")$summary
#' # CTE with growth-weighted temperature average
#' plot(Med_Cc, xlim=c(25, 35))
#' points(x=CTE$TSP.MassWeighted.temperature.mean, y=males/sexed, 
#'          col="red", pch=19)
#' legend("topright", legend=c("CTE with growth-weighted TRN"), 
#'          pch=19, col=c("red"))
#' # CTE with sexualisation TRN and growth-weighted temperature average
#' plot(Med_Cc, xlim=c(25, 35))
#' points(x=CTE$TSP.MassWeighted.STRNWeighted.temperature.mean, y=males/sexed, 
#'          col="red", pch=19)
#' legend("topright", legend=c("CTE with growth-weighted TRN and Sex. TRN"), 
#'        pch=19, col=c("red"))
#' xx <- seq(from=20, to=35, by=0.1)
#' plot(x=xx, 
#'          y=log10(getFromNamespace(".SSM", ns="embryogrowth")(xx, 
#'                             c(fitSTRN$par, fitSTRN$fixed.parameters))[[1]]), 
#'          type="l", bty="n", xlim=c(20, 35), ylim=c(-20, 20), 
#'          xlab="Temperature", ylab="Sexualisation thermal reaction norm (log10)")
#'          
#' # Using only the sexualisation thermal reaction norm within TSP to calculate CTE
#' 
#' Initial_STRN <- resultNest_4p_SSM4p$par[c("DHA", "DHH", "T12H")]
#' Initial_STRN <- structure(c(3678.94960547096, -301.436485427701, 912.595953854977), 
#'                           .Names = c("DHA", "DHH", "T12H"))
#' fp <- c(Rho25=100)
#' fitSTRN_2 <- STRN(Initial_STRN=Initial_STRN, 
#'                 EmbryoGrowthTRN=resultNest_4p_SSM4p, tsd=Med_Cc, 
#'                 Sexed=sexed, Males=males, 
#'                 fixed.parameters=fp,  
#'                 Temperatures="TSP.STRNWeighted.temperature.mean")
#' CTE <- info.nests(NestsResult=resultNest_4p_SSM4p, 
#'                   SexualisationTRN=fitSTRN_2,
#'                   SexualisationTRN.CI="Hessian",
#'                   CI="Hessian", 
#'                   replicate.CI=100, 
#'                   progress=TRUE, 
#'                   warnings=TRUE, 
#'                   out="summary")$summary
#' # CTE with sexualisation TRN
#' plot(Med_Cc, xlim=c(25, 35))
#' points(x=CTE$TSP.STRNWeighted.temperature.mean, y=males/sexed, 
#'          col="red", pch=19)
#' legend("topright", legend=c("CTE with Sexualisation TRN"), 
#'        pch=19, col=c("red"))
#' xx <- seq(from=20, to=35, by=0.1)
#' plot(x=xx, 
#'          y=getFromNamespace(".SSM", ns="embryogrowth")(xx, 
#'                             c(fitSTRN$par, fitSTRN$fixed.parameters))[[1]], 
#'          type="l", bty="n", xlim=c(20, 35), ylim=c(0, 1E-18), 
#'          xlab="Temperature", ylab="Sexualisation thermal reaction norm")

#' }
#' @export

STRN <- function(Initial_STRN=NULL, 
                 fixed.parameters = NULL, 
                 EmbryoGrowthTRN=stop("Embryo Growth Thermal Reaction Norm must be provided"), 
                 tsd=stop("A result from the function tsd() must be provided"),
                 Sexed=NULL, Males=NULL, Females=NULL, 
                 Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean", 
                 SE=TRUE, parallel=TRUE, 
                 itnmax=1000, 
                 method = c("Nelder-Mead","BFGS"), 
                 control=list(trace=1, REPORT=10))
  
{
  
  # Initial_STRN=NULL; fixed.parameters = NULL;EmbryoGrowthTRN=NULL; tsd=NULL; Sexed=NULL; Males=NULL; Females=NULL; Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean"; SE=TRUE; parallel=TRUE; control=list(trace=1, REPORT=10, maxit=1000)
  
  if (!requireNamespace("optimx", quietly = TRUE)) {
    stop("optimx package is absent; Please install it first")
  }
  
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("numDeriv package is absent; Please install it first")
  }
  
#  Initial_STRN=NULL;  EmbryoGrowthTRN=NULL; fixed.parameters = NULL; tsd=NULL;  Sexed=NULL; Males=NULL; Females=NULL;  Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean"; SE=FALSE 
  
  if (is.null(Initial_STRN)) {pSTRN=EmbryoGrowthTRN$par} else {pSTRN=Initial_STRN}
  
  if (is.null(Sexed)) {Sexed <- Males+Females}
  if (is.null(Males)) {Males <- Sexed-Females}
  if (is.null(Females)) {Females <- Sexed-Males}
  
  if (identical(numeric(0), Sexed+Males+Females)) {
    stop("Error in Males, Females or Sexed data")
  }
  
  nm <- names(pSTRN)
  minL <- length(method)
  
  repeat {
    
    L <- list(hessian=FALSE, method=method, 
              par=pSTRN, 
              fixed.parameters=fixed.parameters, 
              fn=getFromNamespace(".STRN_fit", ns="embryogrowth"), 
              EmbryoGrowthTRN=EmbryoGrowthTRN, 
              tsd=tsd, Sexed=Sexed, Males=Males, Temperatures=Temperatures, 
              parallel=parallel, 
              itnmax=itnmax, 
              control=modifyList(list(dowarn=FALSE, follow.on=TRUE, kkt=FALSE), control))
    
    # getFromNamespace(".STRN_fit", ns="embryogrowth")(par=pSTRN, fixed.parameters=fixed.parameters, EmbryoGrowthTRN=EmbryoGrowthTRN, tsd=tsd, Sexed=Sexed, Males=Males, Temperatures=Temperatures)
    
    result <- do.call(getFromNamespace(x="optimx", ns="optimx"), L)
    
    x <- result[minL, nm]
    x <- as.numeric(x)
    names(x) <- nm
    conv <- result[minL, "convcode"]
    value <- result[minL, "value"]
    
    if (conv != 1) break
    pSTRN <- x
    print("Convergence is not achieved. Optimization continues !")
    print(dput(pSTRN))
  }
  
  result <- list(result=result, par = x)
  
  if (SE) {
    mathessian <- try(numDeriv::hessian(getFromNamespace(".STRN_fit", ns="embryogrowth"), 
                                        method="Richardson", 
                                        x=result$par, 
                                        fixed.parameters=fixed.parameters, 
                                        EmbryoGrowthTRN=EmbryoGrowthTRN, 
                                        tsd=tsd, Sexed=Sexed, Males=Males, 
                                        Temperatures=Temperatures
                                        , parallel=parallel), silent=TRUE)
    
    if (inherits(mathessian, "try-error")) {
      res <- rep(NA, length(x))
      names(res) <- names(result$par)
      result$hessian <- NULL
    } else {
      colnames(mathessian) <- rownames(mathessian) <- names(result$par)
      rh <- SEfromHessian(mathessian, hessian=TRUE)
      res <- rh$SE
      result$hessian <- rh$hessian
    }
    
  } else {
    # pas de hessian donc pas de SE
    res<-rep(NA, length(result$par))
    names(res) <- names(result$par)
  }
  
  result$SE <- res 
  
  if (any(is.na(res))) {
    if (all(is.na(res))) {
      print("SE of parameters cannot be estimated.")
      print("Probably the model is badly fitted. Try other initial points.")
    } else {
      print("Probably flat likelihood is observed around some parameters.")
      print("Try using STRN_MHmcmc() function to get the SE of parameters.")
    }
  }
  
  result$data <- list(Sexed=Sexed, Males=Males, Females=Females, 
                      Temperatures="Temperatures", 
                      EmbryoGrowthTRN=EmbryoGrowthTRN, 
                      tsd=tsd)
  
  result$value <- value
  
  k <- length(x)
  n <- length(result$data$Sexed)

  result$AIC <- 2*result$value+2*k
  result$AICc <- result$AIC + (2*k*(k+1))/(n-k-1)
  result$BIC <- 2*result$value + k * log(n)
  
  result$fixed.parameters <- fixed.parameters
  
  class(result) <- "STRN"
  
  return(invisible(result))
}
