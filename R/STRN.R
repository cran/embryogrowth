#' STRN estimates the parameters that best describe the sexualisation thermal reaction norm within the TSP
#' @title Estimate the parameters that best describe the sexualisation thermal reaction norm within the TSP
#' @author Marc Girondot
#' @return The list with object return by optim() 
#' @param Initial_STRN Values for initial model of Sexualisation Thermal Reaction Norm
#' @param EmbryoGrowthTRN The Embryo Growth Thermal Reaction Norm obtained with searchR()
#' @param tsd The model used to predict sex ratio, obtained from tsd()
#' @param Sexed The number of sexed embryos with names identifying timeseries
#' @param Males The number of males embryos with names identifying timeseries
#' @param Females The number of females embryos with names identifying timeseries
#' @param Temperatures The temperature from out of info.nests to be used
#' @param SE Should standard error of parameters be estimated ? TRUE or FALSE
#' @param ... Parameters used for control of optimx()
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
#' # Initial_STRN <- rep(1, 7)
#' # names(Initial_STRN) <- as.character(seq(from=20, to=35, length=7))
#' Initial_STRN <- structure(c(1, 143.248982215757, -25.7029976477549, -0.00489843027318209,
#' -8.94560833594928, 135.781961273868, 71.2176230826628), 
#' .Names = c("20", "22.5", "25", "27.5", "30", "32.5", "35"))
#' males <- c(7, 0, 0, 0, 0, 5, 6, 3, 5, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' names(males) <- rev(rev(names(resultNest_4p$data))[-(1:2)])
#' sexed <- rep(10, length(males))
#' names(sexed) <- rev(rev(names(resultNest_4p$data))[-(1:2)])
#' fitSTRN <- STRN(Initial_STRN, EmbryoGrowthTRN=resultNest_4p, tsd=Med_Cc, 
#' Sexed=sexed, Males=males, 
#' Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean")
#' CTE <- info.nests(NestsResult=resultNest_4p, 
#'  SexualisationTRN=fitSTRN$par, out="summary")$summary
#' plot_add(x=CTE$TSP.MassWeighted.STRNWeighted.temperature.mean, y=males/sexed, 
#'  col="red", pch=19)
#' legend("topright", legend=c("CTE with Sexualisation TRN"), 
#' pch=19, col=c("red"))
#' plotR(parameters=fitSTRN$par, main="Sexualisation TRN")
#' # Initial_STRN <- resultNest_4p$par
#' Initial_STRN <- structure(c(4230.10750319997, 510.543319171189, 1015.78663983953,
#' 118.189709917707), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' males <- c(7, 0, 0, 0, 0, 5, 6, 3, 5, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, NA)
#' names(males) <- rev(rev(names(resultNest_4p$data))[-(1:2)])
#' sexed <- c(rep(10, length(males)-1), NA)
#' names(sexed) <- rev(rev(names(resultNest_4p$data))[-(1:2)])
#' fitSTRN <- STRN(Initial_STRN, EmbryoGrowthTRN=resultNest_4p, tsd=Med_Cc, 
#' Sexed=sexed, Males=males, 
#' Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean")
#' CTE <- info.nests(NestsResult=resultNest_4p, 
#' SexualisationTRN=fitSTRN$par, out="summary")$summary
#' plot(Med_Cc, xlim=c(25, 35))
#' plot_add(x=CTE$TSP.MassWeighted.STRNWeighted.temperature.mean, y=males/sexed, 
#' col="red", pch=19)
#' legend("topright", legend=c("CTE with Sexualisation TRN"), 
#' pch=19, col=c("red"))
#' plotR(parameters=fitSTRN$par, main="Sexualisation TRN")
#' }
#' @export

STRN <- function(Initial_STRN=NULL, 
                 EmbryoGrowthTRN=stop("Embryo Growth Thermal Reaction Norm must be provided"), 
                 tsd=stop("A result from the function tsd() must be provided"),
                 Sexed=NULL, Males=NULL, Females=NULL, 
                 Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean", 
                 SE=FALSE, 
                 ...)
  
{
  
  if (!requireNamespace("optimx", quietly = TRUE)) {
    stop("optimx package is absent; Please install it first")
  }
  
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("numDeriv package is absent; Please install it first")
  }
  
#  Initial_STRN=NULL;  EmbryoGrowthTRN=NULL; tsd=NULL;  Sexed=NULL; Males=NULL; Females=NULL;  Temperatures="TSP.MassWeighted.STRNWeighted.temperature.mean"; SE=FALSE 
  
  if (is.null(Initial_STRN)) {pSTRN=EmbryoGrowthTRN$par} else {pSTRN=Initial_STRN}
  
  if (is.null(Sexed)) {Sexed <- Males+Females}
  if (is.null(Males)) {Males <- Sexed-Females}
  if (is.null(Females)) {Females <- Sexed-Males}
  
  if (identical(numeric(0), Sexed+Males+Females)) {
    stop("Error in Males, Females or Sexed data")
  }
  
  method <- c("Nelder-Mead","BFGS")
  p3p <- list(...)
  # p3p <- list()
  
  repeat {
    
    L <- list(hessian=SE, method=method, 
              par=pSTRN, fn=getFromNamespace(".STRN_fit", ns="embryogrowth"), 
              EmbryoGrowthTRN=EmbryoGrowthTRN, 
              tsd=tsd, Sexed=Sexed, Males=Males, Temperatures=Temperatures, 
              control=modifyList(list(dowarn=FALSE, follow.on=TRUE, kkt=FALSE), p3p))
    
    result <- do.call(getFromNamespace(x="optimx", ns="optimx"), L)
    
    minL <- dim(result)[1]
    nm <- names(pSTRN)
    x <- result[minL, nm]
    x <- as.numeric(x)
    names(x) <- nm
    conv <- result[minL, "convcode"]
    value <- result[minL, "value"]
    
    if (conv==0) break
    pSTRN <- x
    print("Convergence is not achieved. Optimization continues !")
    print(dput(pSTRN))
  }
  
  result <- list(result=result, par = x)
  
  if (SE) {
    mathessian <- try(numDeriv::hessian(getFromNamespace(".STRN_fit", ns="embryogrowth"), 
                                        parameters=result$par, method="Richardson", 
                                        EmbryoGrowthTRN=EmbryoGrowthTRN, 
                                        tsd=tsd, Sexed=Sexed, Males=Males, 
                                        Temperatures=Temperatures), silent=TRUE)
    
    if (inherits(mathessian, "try-error")) {
      res <- rep(NA, length(x))
    } else {
      result$hessian <- mathessian
      inversemathessian <- try(solve(mathessian), silent=TRUE)
      if (inherits(inversemathessian, "try-error")) {
        res <- rep(NA, length(result$par))
      } else {
        res <- diag(inversemathessian)
        res <- ifelse(res<0, NA, sqrt(res))
      }
    }
    
  } else {
    # pas de hessian donc pas de SE
    res<-rep(NA, length(result$par))
  }
  names(res) <- names(result$par)
  result$SE <- res 
  
  if (any(is.na(res))) {
    if (all(is.na(res))) {
      print("SE of parameters cannot be estimated.")
      print("Probably the model is badly fitted. Try other initial points.")
    } else {
      print("Probably flat likelihood is observed around some parameters.")
      print("Try using GRTRN_MHmcmc() function to get the SE of parameters.")
    }
  }
  
  result$value <- value
  result$AIC <- 2 * value + 2* length(x)
  result$data <- list(Sexed=Sexed, Males=Males, Females=Females, 
                      Temperatures="Temperatures", 
                      EmbryoGrowthTRN=EmbryoGrowthTRN, 
                      tsd=tsd)
  
  class(result) <- "STRN"
  
  return(invisible(result))
}
