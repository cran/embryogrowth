#' STRN estimates the parameters that best describe the sexualisation thermal reaction norm within the TSP
#' @title Estimate the parameters that best describe the sexualisation thermal reaction norm within the TSP
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The list with object return by optim() 
#' @param Initial_STRN Values for initial model of Sexualisation Thermal Reaction Norm or tsd model
#' @param fixed.parameters Value for Sexualisation Thermal Reaction Norm or tsd model that will not be changed
#' @param EmbryoGrowthTRN The Embryo Growth Thermal Reaction Norm obtained with searchR()
#' @param embryo.stages The embryo stages. At least TSP.borders stages must be provided to estimate TSP borders. See note.
#' @param TSP.borders The limits of TSP in stages. See embryo.stages parameter.
#' @param TSP.begin Where TSP begin during the stage of beginning? In relative proportion of the stage.
#' @param TSP.end Where TSP begin during the stage of ending? In relative proportion of the stage.
#' @param tsd The model used to predict sex ratio, obtained from tsd()
#' @param equation If tsd parameter is not provided, equation and parameters in Initial_STRN for tsd model must be provided.
#' @param Sexed The number of sexed embryos with names identifying timeseries
#' @param Males The number of males embryos with names identifying timeseries
#' @param Females The number of females embryos with names identifying timeseries
#' @param sexratio The sex ratio to be used
#' @param method Methods to be used with optimx
#' @param fill See info.nests()
#' @param itnmax Maximum number of iterations for each method; if 0, just return the likelihood
#' @param parallel Should parallel computing for info.nests() be used
#' @param control List for control parameters for optimx
#' @param zero The value to replace a null sex ratio
#' @param verbose If TRUE, will show all intermediate parameters during fit
#' @param hessian If TRUE, the Hessian approximation is estimated atthe end of the fit.
#' @description Estimate the parameters that best describe the sexualisation thermal reaction norm within the TSP.\cr
#' The sexratio parameter is a character string which can be:\cr
#' \itemize{
#'   \item \code{TSP.TimeWeighted.sexratio.mean} Sex ratio based on average temperature during the TSP
#'   \item \code{TSP.GrowthWeighted.sexratio.mean} Sex ratio based on average temperature weighted by the actual growth during the TSP
#'   \item \code{TSP.TimeWeighted.GrowthRateWeighted.sexratio.mean} Sex ratio based on average temperature weighted by the growth rate during the TSP
#'   \item \code{TSP.TimeWeighted.STRNWeighted.sexratio.mean} Sex ratio based on average temperature weighted by the thermal reaction norm of sexualization during the TSP
#'   \item \code{TSP.GrowthWeighted.STRNWeighted.sexratio.mean} Sex ratio based on average temperature weighted by the actual growth and thermal reaction norm of sexualization during the TSP
#'   \item \code{TSP.TimeWeighted.GrowthRateWeighted.STRNWeighted.sexratio.mean} Sex ratio based on average temperature weighted by the growth rate and the thermal reaction norm of sexualization during the TSP
#'   \item \code{MiddleThird.TimeWeighted.sexratio.mean} Sex ratio based on average temperature during the middle third incubation
#'   \item \code{MiddleThird.GrowthWeighted.sexratio.mean} Sex ratio based on average temperature weighted by actual growth during the middle third incubation
#'   \item \code{MiddleThird.TimeWeighted.GrowthRateWeighted.sexratio.mean} Sex ratio based on average temperature weighted by growth rate during the middle third incubation
#'   \item \code{TimeWeighted.sexratio.mean} Sex ratio based on average temperature during all incubation
#'   \item \code{GrowthWeighted.sexratio.mean} Sex ratio based on average temperature weighted by actual growth during all incubation
#'   \item \code{TimeWeighted.GrowthRateWeighted.sexratio.mean} Sex ratio based on average temperature weighted by growth rate during all incubation
#'   \item \code{TSP.PM.TimeWeighted.mean} Average sex ratio based on temperature during the TSP
#'   \item \code{TSP.PM.GrowthWeighted.mean} Average sex ratio based on temperature weighted by the actual growth during the TSP
#'   \item \code{TSP.PM.TimeWeighted.GrowthRateWeighted.mean} Average sex ratio based on temperature weighted by the growth rate during the TSP
#'   }
#' If information for sex is not known for some timeseries, set NA for Sexed.\cr
#' Sexed, Males and Females must be vectors with names. The names must be the same as 
#' the names of timeseries of temperatures in EmbryoGrowthTRN.\cr
#' Only two of these 3 parameters are required: Males, Females and Sexed\cr
#' Note: four species have predefined embryo stages. embryo.stages parameter can take the values:\cr
#' \itemize{
#'   \item \code{Caretta caretta.SCL}
#'   \item \code{Chelonia mydas.SCL}
#'   \item \code{Emys orbicularis.SCL}
#'   \item \code{Emys orbicularis.mass}
#'   \item \code{Podocnemis expansa.SCL}
#'   \item \code{Lepidochelys olivacea.SCL}
#'   \item \code{Generic.ProportionDevelopment}
#' }
#' A fifth name \code{fitted} must be used when limits of TSP are fitted using \code{BeginTSP} and \code{EndTSP} parameters.\cr
#' The parameters that can be used in STRN are:\cr
#' \code{BeginTSP}, \code{EndTSP} are the logit of the proportion of development; \cr
#' To ensure that \code{BeginTSP} < \code{EndTSP}, it is better to use:\cr
#' \code{BeginTSP}, \code{LengthTSP} and then \code{EndTSP} is estimated using \code{BeginTSP} + \code{abs(LengthTSP)}\cr
#' \code{DHA}, \code{DHH}, \code{T12H} are the SSM parameters of sexualisation thermal reaction norm; \cr
#' \code{dbeta_mu}, \code{dbeta_v} are the beta mean and variance of the impact of sexualisation according to TSP progress.\cr
#' Or any parameter that can be used in a TSD model.
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' MedIncubation_Cc <- subset(DatabaseTSD, Species=="Caretta caretta" & 
#' RMU=="Mediterranean" & Sexed!=0)
#' Med_Cc <- tsd(males=MedIncubation_Cc$Males, 
#'              females=MedIncubation_Cc$Females, 
#'              temperatures=MedIncubation_Cc$Incubation.temperature, 
#'              par=c(P=29.5, S=-0.1))
#' plot(Med_Cc, xlim=c(25, 35))
#' males <- c(7, 0, 0, 0, 0, 5, 6, 3, 5, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' names(males) <- rev(rev(names(resultNest_4p_SSM$data))[-(1:2)])
#' sexed <- rep(10, length(males))
#' names(sexed) <- rev(rev(names(resultNest_4p_SSM$data))[-(1:2)])
#' 
#' Initial_STRN <- c('DHA' = 1174.6461503413307, 
#'                   'DHH' = 2001.0619192107047, 
#'                   'T12H' = 3731.353104743393)
#' fp <- c(Rho25=100)
#' fitSTRN <- STRN(Initial_STRN=Initial_STRN, 
#'                 EmbryoGrowthTRN=resultNest_4p_SSM, 
#'                 tsd=Med_Cc, 
#'                 embryo.stages="Caretta caretta.SCL", 
#'                 Sexed=sexed, Males=males, 
#'                 fixed.parameters=fp, 
#'                 sexratio="TSP.GrowthWeighted.STRNWeighted.sexratio.mean")
#' plotR(fitSTRN, curve ="ML", ylim=c(0,2))
#' plotR(fitSTRN)
#' out <- info.nests(NestsResult=resultNest_4p_SSM, 
#'                   SexualisationTRN=fitSTRN,
#'                   SexualisationTRN.CI="Hessian",
#'                   embryo.stages="Caretta caretta.SCL", 
#'                   GTRN.CI="Hessian", 
#'                   tsd=Med_Cc, 
#'                   tsd.CI="Hessian", 
#'                   replicate.CI=100, 
#'                   progressbar=TRUE, 
#'                   warnings=TRUE, 
#'                   out="summary")$summary
#' # CTE with growth-weighted temperature average
#' plot(Med_Cc, xlim=c(25, 35))
#' points(x=out[, "TSP.GrowthWeighted.STRNWeighted.temperature.mean"], y=males/sexed, 
#'          col="red", pch=19)
#' legend("topright", legend=c("CTE with growth-weighted and Sexualization TRN"), 
#'          pch=19, col=c("red"))
#'        
#' #  Fit the beginning and end of TSP
#' 
#' Initial_STRN <- c('BeginTSP' = invlogit(0.33), 
#'                   'EndTSP' = invlogit(0.66))
#' fp <- NULL
#' fitSTRN <- STRN(Initial_STRN=Initial_STRN, 
#'                 EmbryoGrowthTRN=resultNest_4p_SSM, 
#'                 tsd=Med_Cc, 
#'                 embryo.stages="fitted", 
#'                 Sexed=sexed, Males=males, 
#'                 fixed.parameters=fp, 
#'                 sexratio="TSP.TimeWeighted.GrowthRateWeighted.STRNWeighted.sexratio.mean")
#' invlogit(fitSTRN$par)
#' invlogit(fitSTRN$par-2*fitSTRN$SE)
#' invlogit(fitSTRN$par+2*fitSTRN$SE)
#'
#' Initial_STRN <- c('dbeta_mu' = logit(0.5), 
#'                   'dbeta_v' = 1/12)
#' fp <- NULL
#' fitSTRN <- STRN(Initial_STRN=Initial_STRN, 
#'                 EmbryoGrowthTRN=resultNest_4p_SSM, 
#'                 tsd=Med_Cc, 
#'                 embryo.stages="Caretta caretta.SCL", 
#'                 Sexed=sexed, Males=males, 
#'                 fixed.parameters=fp, 
#'                 sexratio="TSP.TimeWeighted.GrowthRateWeighted.STRNWeighted.sexratio.mean")
#'                 mu <- invlogit(fitSTRN$par["dbeta_mu"]), 
#'                 v <- abs(fitSTRN$par["dbeta_v"])
#'                 shape1 <- mu * (((mu * (1 - mu))/v) - 1)
#'                 shape2 <- shape1 * (1 - mu)/mu
#' plot(seq(from=0, to=1, length.out=100), 
#'      dbeta(seq(from=0, to=1, length.out=100), 
#'                shape1=shape1, shape2=shape2), 
#'      type="l", xlab="Progress of TSP", 
#'      ylab="Force of sexualisation", bty="n", ylim=c(0, 6), las=1)
#'                
#'  Initial_STRN <- c('dbeta_mu' = 7.2194053298953236, 
#'                    'dbeta_v' = 0.00050390986089928467)
#' fp <- NULL
#' fitSTRN <- STRN(Initial_STRN=Initial_STRN                                                 , 
#'                 EmbryoGrowthTRN=resultNest_4p_SSM                                         , 
#'                 tsd=Med_Cc                                                                , 
#'                 embryo.stages="Caretta caretta.SCL"                                       , 
#'                 Sexed=sexed                                                               , 
#'                 Males=males                                                               , 
#'                 fixed.parameters=fp                                                       , 
#'                 sexratio="TSP.TimeWeighted.GrowthRateWeighted.STRNWeighted.sexratio.mean" )
#'                 
#'                 mu <- invlogit(fitSTRN$par["dbeta_mu"]), 
#'                 v <- abs(fitSTRN$par["dbeta_v"])
#'                 shape1 <- mu * (((mu * (1 - mu))/v) - 1)
#'                 shape2 <- shape1 * (1 - mu)/mu
#'                 
#' plot(seq(from=0, to=1, length.out=100), 
#'      dbeta(seq(from=0, to=1, length.out=100), 
#'                shape1=shape1, shape2=shape2), 
#'      type="l", xlab="Progress of TSP", 
#'      ylab="Force of sexualisation", bty="n", ylim=c(0, 0.04), las=1)
#'                
#' Initial_STRN <- c('dbeta_mu' = logit(0.5), 
#'                   'dbeta_v' = 1/12)
#' L <- STRN(Initial_STRN=NULL                                                        , 
#'           fixed.parameters=Initial_STRN                                            , 
#'           EmbryoGrowthTRN=resultNest_4p_SSM                                        , 
#'           tsd=Med_Cc                                                               , 
#'           embryo.stages="Caretta caretta.SCL"                                                   , 
#'           Sexed=sexed                                                              ,
#'           Males=males                                                              , 
#'           sexratio="TSP.TimeWeighted.GrowthRateWeighted.STRNWeighted.sexratio.mean")
#' Initial_STRN <- c('dbeta_mu' = logit(0.6), 
#'                   'dbeta_v' = 1/12)
#' L <- STRN(Initial_STRN=NULL                                                        , 
#'           fixed.parameters=Initial_STRN                                            , 
#'           EmbryoGrowthTRN=resultNest_4p_SSM                                        , 
#'           tsd=Med_Cc                                                               , 
#'           embryo.stages="Caretta caretta.SCL"                                                   , 
#'           Sexed=sexed                                                              ,
#'           Males=males                                                              , 
#'           sexratio="TSP.TimeWeighted.GrowthRateWeighted.STRNWeighted.sexratio.mean")
#'  Initial_STRN <- c('dbeta_mu' = 7.2192972077000004, 
#'                     'dbeta_v' = 0.00050396969999999997)
#' L <- STRN(Initial_STRN=NULL                                                        , 
#'           fixed.parameters=Initial_STRN                                            , 
#'           EmbryoGrowthTRN=resultNest_4p_SSM                                        , 
#'           tsd=Med_Cc                                                               , 
#'           embryo.stages="Caretta caretta.SCL"                                      , 
#'           Sexed=sexed                                                              ,
#'           Males=males                                                              , 
#'           sexratio="TSP.TimeWeighted.GrowthRateWeighted.STRNWeighted.sexratio.mean")
#'                 mu <- invlogit(fitSTRN$par["dbeta_mu"]), 
#'                 v <- abs(fitSTRN$par["dbeta_v"])
#'                 shape1 <- mu * (((mu * (1 - mu))/v) - 1)
#'                 shape2 <- shape1 * (1 - mu)/mu
#'  
#'    tsp_progress <- seq(from=0, to=1, length.out=100)     
#'   plot(tsp_progress, 
#'      dbeta(tsp_progress, 
#'            shape1=shape1, shape2=shape2), 
#'        type="l", xlab="Progress of TSP", 
#'        ylab="Force of sexualisation", bty="n", ylim=c(0, 0.04), las=1)
#'  segments(x0=0, x1=1, y0=0, y1=0, lty=2)
#' }
#' @export

STRN <- function(EmbryoGrowthTRN=stop("Embryo Growth Thermal Reaction Norm must be provided"),
                 Initial_STRN=NULL                                                           , 
                 fixed.parameters = NULL                                                     , 
                 TSP.borders=NULL                                                            , 
                 embryo.stages=NULL                                                          , 
                 TSP.begin=0                                                                 , 
                 TSP.end=0.5                                                                 , 
                 tsd=NULL                                                                    ,
                 equation="logistic"                                                         , 
                 Sexed=NULL                                                                  , 
                 Males=NULL                                                                  , 
                 Females=NULL                                                                , 
                 sexratio="TSP.TimeWeighted.GrowthRateWeighted.STRNWeighted.sexratio.mean"   , 
                 fill = 60                                                                   ,
                 parallel=TRUE                                                               , 
                 itnmax=1000                                                                 , 
                 method = c("Nelder-Mead", "BFGS")                                           , 
                 control=list(trace=1, REPORT=10)                                            , 
                 zero=1E-9                                                                   , 
                 verbose=FALSE                                                               , 
                 hessian=TRUE                                                                ) {
  
  # Initial_STRN=NULL; fixed.parameters = NULL; EmbryoGrowthTRN=NULL; TSP.borders=NULL; TSP.begin=0; TSP.end=0.5; embryo.stages=NULL; tsd=NULL; equation="logistic"; Sexed=NULL; Males=NULL; Females=NULL; sexratio="TSP.GrowthWeighted.STRNWeighted.sexratio"; SE=TRUE; parallel=TRUE; itnmax=1000; method = c("Nelder-Mead","BFGS"); control=list(trace=1, REPORT=10); zero=1E-9
  
  if (is.null(embryo.stages)) {
    stop("embryo.stages must be defined")
  }
  
  if (!requireNamespace("optimx", quietly = TRUE)) {
    stop("optimx package is absent; Please install it first")
  }
  
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("numDeriv package is absent; Please install it first")
  }
  
  pSTRN <- Initial_STRN
  
  if (is.null(Sexed)) {Sexed <- Males+Females}
  if (is.null(Males)) {Males <- Sexed-Females}
  if (is.null(Females)) {Females <- Sexed-Males}
  
  if (identical(numeric(0), Sexed+Males+Females)) {
    stop("Error in Males, Females or Sexed data")
  }
  
  nm <- names(pSTRN)
  minL <- length(method)
  
  if ((is.null(pSTRN)) | (itnmax == 0)) {
    x <- NULL
    value <- getFromNamespace(".STRN_fit", ns="embryogrowth")(par=pSTRN                          , 
                                                              fixed.parameters=fixed.parameters  ,
                                                              equation=equation                  , 
                                                              tsd=tsd                            , 
                                                              Sexed=Sexed                        , 
                                                              Males=Males                        , 
                                                              sexratio=sexratio                  , 
                                                              parallel=parallel                  , 
                                                              zero=zero                          , 
                                                              NestsResult = EmbryoGrowthTRN      , 
                                                              embryo.stages=embryo.stages        , 
                                                              TSP.borders=TSP.borders            , 
                                                              TSP.begin=TSP.begin                , 
                                                              TSP.end=TSP.end                    , 
                                                              fill=fill                          ,
                                                              out="likelihood"                   ,
                                                              verbose=verbose                    )
    result <- list(value=value)
    result$par <- pSTRN
  } else {
    
    repeat {
      
      L <- list(hessian=hessian                                                            , 
                method=method                                                              , 
                par=pSTRN                                                                  , 
                fixed.parameters=fixed.parameters                                          , 
                fn=getFromNamespace(".STRN_fit", ns="embryogrowth")                        , 
                equation=equation                                                          , 
                tsd=tsd                                                                    , 
                Sexed=Sexed                                                                , 
                Males=Males                                                                , 
                sexratio=sexratio                                                          , 
                parallel=parallel                                                          , 
                zero=zero                                                                  ,
                verbose=verbose                                                            , 
                NestsResult = EmbryoGrowthTRN                                              , 
                embryo.stages=embryo.stages                                                , 
                TSP.borders=TSP.borders                                                    , 
                TSP.begin=TSP.begin                                                        , 
                TSP.end=TSP.end                                                            ,
                fill=fill                                                                  ,
                itnmax=itnmax                                                              , 
                control=modifyList(list(dowarn=FALSE, follow.on=TRUE, kkt=FALSE), control) )
      
      
      result <- do.call(getFromNamespace(x="optimx", ns="optimx"), L)
      
      colnames(result) <- c(nm, colnames(result)[(length(nm)+1): ncol(result)])
      
      x <- unlist(result[minL, nm, drop=TRUE])
      # x <- as.numeric(x)
      # names(x) <- nm
      conv <- result[minL, "convcode"]
      value <- result[minL, "value"]
      
      if (conv != 1) break
      pSTRN <- x
      print("Convergence is not achieved. Optimization continues !")
      print(dput(pSTRN))
    }
    
    result <- list(result=result, par = x, hessian=attr(result, "details")[nrow(result), ]$nhatend)
    result$value <- value
    k <- length(x)
    
    if (!hessian) {
      result$hessian <- matrix(NA, ncol=k, nrow=k)
    }
    
    colnames(result$hessian) <- rownames(result$hessian) <- names(x)
    
    if (all(!is.na(result$hessian))) {
      result$SE <- SEfromHessian(result$hessian, hessian=FALSE)
    } else {
      result$SE <- rep(NA, length(nm))
      names(result$SE) <- nm
    }
    
    if (any(is.na(result$SE)) & hessian) {
      print("SE of parameters cannot be estimated.")
      print("Probably the model is badly fitted. Try other initial points or MCMC.")
    }
  }
  
  k <- length(x)
  
  #22/2/2020, ajout de TSP; pas sûr si ca sert
  result$Sexed <- Sexed
  result$Males <- Males
  result$Females <- Females
  result$tsd <- tsd
  result$equation <- equation
  
  result$EmbryoGrowthTRN <- EmbryoGrowthTRN
  result$embryo.stages <- embryo.stages                       
  result$TSP.borders <- TSP.borders                           
  result$TSP.begin <- TSP.begin                         
  result$TSP.end <- TSP.end
  result$fill <- fill
  result$zero <- zero
  
  n <- sum(result$Sexed, na.rm = TRUE)
  
  result$AIC <- 2*result$value+2*k
  result$AICc <- result$AIC + (2*k*(k+1))/(n-k-1)
  result$BIC <- 2*result$value + k * log(n)
  
  result$fixed.parameters <- fixed.parameters
  
  result$sexratio <- sexratio
  
  
  
  result <- addS3Class(result, "STRN")
  
  return(invisible(result))
}
