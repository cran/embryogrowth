#' likelihoodR estimates the likelihood of a set of parameters for nest incubation data 
#' @title Estimate the likelihood of a set of parameters for nest incubation data
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A result object
#' @param result A object obtained after searchR or likelihoodR
#' @param parameters A set of parameters
#' @param fixed.parameters A set of parameters that will not be changed
#' @param temperatures Timeseries of temperatures
#' @param integral Function used to fit embryo growth: integral.Gompertz, integral.exponential or integral.linear
#' @param derivate Function used to fit embryo growth: dydt.Gompertz, dydt.exponential or dydt.linear. It will replace the one in NestsResult.
#' @param hatchling.metric Mean and SD of size of hatchlings
#' @param M0 Measure of hatchling size or mass proxi at laying date
#' @param weight A named vector of the weight for each nest for likelihood estimation
#' @param hessian If TRUE, the hessian matrix is estimated and the SE of parameters estimated.
#' @param parallel If true, try to use several cores using parallel computing.
#' @param echo If FALSE, does not display the result.
#' @description Estimate the likelihood of a set of parameters for nest incubation data
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
#' # K or rK are not used for integral.linear or integral.exponential
#' LresultNest_4p <- likelihoodR(parameters=x, fixed.parameters=pfixed,  
#' 	temperatures=formated, integral=integral.Gompertz, M0=1.7,  
#' 	hatchling.metric=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p_SSM)
#' LresultNest_4p <- likelihoodR(result=resultNest_4p_SSM)
#' }
#' @export


likelihoodR <-
  function(result=NULL               , 
           parameters=NULL           , 
           fixed.parameters=NULL     , 
           temperatures=NULL         , 
           integral=NULL             , 
           derivate=NULL             , 
           hatchling.metric=NULL     , 
           M0=NULL                   , 
           hessian=FALSE             , 
           weight= NULL              , 
           parallel=TRUE             ,
           echo=TRUE                 ) {
    
    # result=NULL; parameters=NULLfixed.parameters=NULL; temperatures=NULL; echo=TRUE;integral=NULL; derivate=NULL; hatchling.metric=NULL; M0=NULL; hessian=FALSE; weight= NULL; parallel=TRUE
    # result=resultNest_4p
    # parameters=structure(c(0.840429263837085, 1.15929190273597, 1.57697305167496, 2.06657337881587, 2.43555105967532, 2.24337041114994, 1.44633577322632), .Names = c("296", "298", "300", "302", "304", "306", "308"))
    
    
    
    if (!is.null(result)) {
      if (is.null(temperatures)) {
        temperatures <- result[["data"]]
        if (inherits(temperatures, "Nests")) {
          temperatures <- FormatNests(temperatures)
        }
      }
      if (is.null(integral)) integral <- result[["integral"]]
      if (is.null(derivate)) derivate <- result[["derivate"]]
      if (is.null(hatchling.metric)) hatchling.metric <- result[["hatchling.metric"]]
      if (is.null(M0)) M0 <- result[["M0"]]
      if (is.null(fixed.parameters)) fixed.parameters <- result[["fixed.parameters"]]
      if (is.null(parameters)) parameters <- result[["par"]]
      # correction le 29/1/2014
      if (is.null(weight)) weight <- result[["weight"]]
      if (is.null(hatchling.metric)) hatchling.metric <- result$hatchling.metric
    } 
    if (is.numeric(hatchling.metric)) {
      hatchling.metric <- data.frame(Mean=rep(hatchling.metric["Mean"], temperatures[["IndiceT"]]["NbTS"]), SD=rep(hatchling.metric["SD"], temperatures[["IndiceT"]]["NbTS"]), row.names=temperatures$Names)
    }
    
    NbTS <- temperatures$IndiceT["NbTS"]
    
    if (is.null(weight)) {
      par <- setNames(rep(1, NbTS), temperatures$Names)
    } else {
      
      if (any(is.na(weight))) {
        par <- setNames(rep(1, NbTS), temperatures$Names)
      } else {
        
        if (is.list(weight)) weight <- weight$weight
        
        if (length(setdiff(names(weight), temperatures$Names))==0) {
          par <- weight
        } else {
          warning("Check the weights")
          return(invisible())
        }
      }
    }
    
    weight <- par
    
    # te_st si tous sont l
    if (length(setdiff(temperatures$Names, names(weight)))!=0) {
      warning("The weight parameter must define weight for each nest.")
      warning(paste("check", setdiff(temperatures$Names, names(weight)), "nests"))
      return(invisible())	
    }	
    
    
    # Un paramtre ne peut pas tre indique en fixe et en fite - 22/7/2012	
    
    if (length(intersect(names(parameters), names(fixed.parameters)))!=0) {
      warning("A parameter cannot be fixed and fitted at the same time !")
      return(invisible())
    }
    
    
    for (j in temperatures$Names) temperatures$Nests[[j]]$data[1, "Mass"] <- M0
    
    resultnest <- list(par=parameters)
    
    resultnest$value <- info.nests(parameters=parameters, temperatures=temperatures, integral=integral, derivate=derivate, 
                                   hatchling.metric=hatchling.metric, M0=M0, fixed.parameters=fixed.parameters, weight=weight)
    
    if(echo) {
      
      print(paste("-Ln L=",resultnest$value))
      
    }
    
    if (hessian) {
      
      print("Estimate the SE for parameters at that point")
      
      mathessian <- try(hessian(info.nests, parameters=result$par, method="Richardson", temperatures=temperatures, integral=integral, derivate=derivate, 
                                hatchling.metric=hatchling.metric, M0=M0, fixed.parameters=fixed.parameters, weight=weight), silent=TRUE)
      
      if (inherits(mathessian, "try-error")) {
        res<-rep(NA, length(parameters))
        resultnest$hessian <- NA
        print("SE of parameters cannot be estimated.")
        print("Probably the model is badly fitted. Try using searchR() before.")
        
      } else {
        
        
        resultnest$hessian <- mathessian
        
        inversemathessian <- try(solve(mathessian), silent=TRUE)
        
        if (inherits(inversemathessian, "try-error")) {
          res <- -1
        } else {
          res <- diag(inversemathessian)
        }
        
        
        # Je gre plus correctement les erreurs - 17/7/2012
        
        neg=any(res<0)
        if (!neg) {
          res=sqrt(res)
        } else {
          res<-rep(NA, length(parameters))
          print("SE of parameters cannot be estimated.")
          if (inherits(inversemathessian, "try-error")) {
            print("Probably the model is badly fitted. Try using searchR() before.")
          } else {
            print("Probably flat likelihood is observed around some parameters.")
            print("Try using MHmcmc() function to get the SE of parameters.")
          }
        }
        # fin du te_st sur la premire erreur
      }
      
    } else {
      # pas de hessian donc pas de SE
      res<-rep(NA, length(parameters))
    }
    
    
    names(res) <- names(parameters)
    
    resultnest$SE <- res 
    
    resultnest$data <- temperatures
    k <- length(parameters)
    n <- unname(resultnest$data$IndiceT["NbTS"])
    # Avant *5. Correction du 17/7/2012
    if ((all(!is.na(hatchling.metric[, "SD"]))) | (!is.na(resultnest$par["SD"]))) {
    resultnest$AIC <- 2*resultnest$value+k
    resultnest$AICc <- result$AIC + (2*k*(k+1))/(n-k-1)
    resultnest$BIC <- 2*result$value + k * log(n)
    } else {
      resultnest$AIC <- n*log(resultnest$value/n)+2*k
      resultnest$AICc <- resultnest$AIC + (2*k*(k+1))/(n-k-1)
      resultnest$BIC <- n*log(resultnest$value/n) + k * log(n)
    }
    
    
    resultnest$hatchling.metric <- hatchling.metric
    
    resultnest$integral <- integral
    resultnest$derivate <- derivate
    
    resultnest$M0 <- M0
    
    # 29/1/2014
    resultnest$weight <- weight
    
    # Je stocke aussi les paramtres fixe-16/7/2012
    resultnest$fixed.parameters <- fixed.parameters
    
    resultnest <- addS3Class(resultnest, "NestsResult")
    
    return(resultnest)
    
  }
