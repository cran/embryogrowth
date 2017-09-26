#' tsd estimates the parameters that best describe temperature-dependent sex determination
#' @title Estimate the parameters that best describe temperature-dependent sex determination
#' @author Marc Girondot
#' @return A list the pivotal temperature, transitional range of temperatures and their SE
#' @param males A vector with male numbers
#' @param females A vector with female numbers
#' @param N A vector with total numbers
#' @param temperatures The constant incubation temperatures used to fit sex ratio
#' @param durations The duration of incubation or TSP used to fit sex ratio
#' @param df A dataframe with at least two columns named males, females or N and temperatures, Incubation.temperature or durations column
#' @param l Sex ratio limits to define TRT are l and 1-l (see Girondot, 1999)
#' @param parameters.initial Initial values for P, S or K search as a vector, ex. c(P=29, S=-0.3)
#' @param fixed.parameters Parameters that will not be changed
#' @param males.freq If TRUE data are shown as males frequency
#' @param equation Could be "logistic", "Hill", "Richards", "Hulin", "Double-Richards" or "GSD"
#' @param replicate.CI Number of replicates to estimate confidence intervals
#' @param range.CI The range of confidence interval for estimation, default=0.95
#' @param print Do the results must be printed at screen? TRUE (default) or FALSE
#' @description Estimate the parameters that best describe the thermal reaction norm for sex ratio when temperature-dependent sex determination occurs.\cr
#' It can be used also to evaluate the relationship between incubation duration and sex ratio.\cr
#' The parameter l was defined in Girondot (1999). The TRT is defined from the difference between the two boundary temperatures giving sex ratios of \eqn{l} and \eqn{1 - l}, respectively:\cr
#' \deqn{TRT_{l}=abs\left ( S\: K_{l} \right )}{TRTl = abs( S.Kl )}
#' where \eqn{K_{l}}{Kl} is a constant equal to \eqn{2\: log\left ( \frac{l}{1-l} \right )}{2 log(l / ( 1 - l))}.\cr
#' In Girondot (1999), l was 0.05 and then the TRT was defined as being the range of temperatures producing from 5\% to 95\% of each sex.\cr
#' Models Richards, Double-Richards and Hulin are particularly sensitive for K parameters. If you want estimate 
#' confidence interval for TRT using these models, it is better to use mcmc.
#' @references Girondot, M. 1999. Statistical description of temperature-dependent sex determination using maximum likelihood. Evolutionary Ecology Research, 1, 479-486.
#' @references Godfrey, M.H., Delmas, V., Girondot, M., 2003. Assessment of patterns of temperature-dependent sex determination using maximum likelihood model selection. Ecoscience 10, 265-272.
#' @references Hulin, V., Delmas, V., Girondot, M., Godfrey, M.H., Guillon, J.-M., 2009. Temperature-dependent sex determination and global change: are some species at greater risk? Oecologia 160, 493-506.
#' @family Functions for temperature-dependent sex determination
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' CC_AtlanticSW <- subset(DatabaseTSD, RMU=="Atlantic, SW" & 
#'                           Species=="Caretta caretta" & (!is.na(Sexed) & Sexed!=0))
#' tsdL <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="logistic", replicate.CI=NULL))
#' tsdH <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Hill", replicate.CI=NULL))
#' tsdR <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Richards", replicate.CI=NULL))
#' tsdDR <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Double-Richards", replicate.CI=NULL))
#' gsd <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="GSD", replicate.CI=NULL))
#' compare_AIC(Logistic_Model=tsdL, Hill_model=tsdH, Richards_model=tsdR, 
#'                DoubleRichards_model=tsdDR, GSD_model=gsd)
#' compare_AICc(Logistic_Model=tsdL, Hill_model=tsdH, Richards_model=tsdR, 
#'                DoubleRichards_model=tsdDR, GSD_model=gsd, factor.value = -1)
#' compare_BIC(Logistic_Model=tsdL, Hill_model=tsdH, Richards_model=tsdR, 
#'                DoubleRichards_model=tsdDR, GSD_model=gsd, factor.value = -1)
#' ##############
#' eo <- subset(DatabaseTSD, Species=="Emys orbicularis", c("Males", "Females", 
#'                                        "Incubation.temperature"))
#'                                        
#' eo_Hill <- with(eo, tsd(males=Males, females=Females, 
#'                                        temperatures=Incubation.temperature,
#'                                        equation="Hill"))
#' eo_Hill <- tsd(df=eo, equation="Hill", replicate.CI=NULL)
#' eo_logistic <- tsd(eo, replicate.CI=NULL)
#' eo_Richards <- with(eo, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature, 
#'                                  equation="Richards", replicate.CI=NULL))
#' ### The Hulin model is a modification of Richards (See Hulin et al. 2009)
#' par <- eo_Richards$par
#' names(par)[which(names(par)=="K")] <- "K2"
#' par <- c(par, K1=0)
#' eo_Hulin <- with(eo, tsd(males=Males, females=Females, 
#'                                  parameters.initial=par, 
#'                                  temperatures=Incubation.temperature, 
#'                                  equation="Hulin", replicate.CI=NULL))
#' ### The Double-Richards model is a Richards model with K1 and K2 using the two values
#' ### below and above P
#' par <- eo_Richards$par
#' names(par)[which(names(par)=="K")] <- "K2"
#' par <- c(par, K1=as.numeric(par["K2"])-0.1)
#' par["K1"] <- par["K1"]-0.1
#' eo_Double_Richards <- with(eo, tsd(males=Males, females=Females,
#'                                  parameters.initial=par,
#'                                  temperatures=Incubation.temperature,
#'                                  equation="Double-Richards", replicate.CI=NULL))
#' compare_AIC(Logistic=eo_logistic, Hill=eo_Hill, Richards=eo_Richards, 
#'              Hulin=eo_Hulin, Double_Richards=eo_Double_Richards)
#' ### Note the asymmetry of the Double-Richards model
#' predict(eo_Double_Richards, 
#'        temperatures=c(eo_Double_Richards$par["P"]-0.2, eo_Double_Richards$par["P"]+0.2))
#' predict(eo_Double_Richards)
#' ### It can be used also for incubation duration
#' CC_AtlanticSW <- subset(DatabaseTSD, RMU=="Atlantic, SW" & 
#'                           Species=="Caretta caretta" & Sexed!=0)
#' tsdL_IP <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  durations=IP.mean, 
#'                                  equation="logistic", replicate.CI=NULL))
#' plot(tsdL_IP, xlab="Incubation durations in days")
#' }
#' @export


tsd <- function(df=NULL, males=NULL, females=NULL, N=NULL, 
                temperatures=NULL, 
                durations=NULL,
                l=0.05, parameters.initial=c(P=NA, S=-2, K=0, K1=1, K2=0), 
                males.freq=TRUE, 
                fixed.parameters=NULL, 
                equation="logistic", replicate.CI=10000, range.CI=0.95, 
                print=TRUE) {
  
  # df=NULL; males=NULL; females=NULL; N=NULL; temperatures=NULL; durations=NULL; l=0.05; parameters.initial=c(P=NA, S=-0.5, K=0, K1=1, K2=0); males.freq=TRUE; fixed.parameters=NULL; equation="logistic"; replicate.CI=10000; range.CI=0.95; print=TRUE
  # CC_AtlanticSW <- subset(DatabaseTSD, RMU=="Atlantic, SW" & Species=="Caretta caretta" & Sexed!=0)
  # males=CC_AtlanticSW$Males; females=CC_AtlanticSW$Females; temperatures=CC_AtlanticSW$Incubation.temperature-CC_AtlanticSW$Correction.factor; equation="Richards"
  
  equation <- tolower(equation)
  
  if (!is.null(df)) {
    if (class(df)!="data.frame" & class(df)!="matrix") {
      stop("df parameter must be a data.frame or a matrix")
    }
    
    namesdf <- tolower(colnames(df))
    males <- NULL
    females <- NULL
    N <- NULL
    
    if (any(namesdf=="males")) males <- df[, which(namesdf=="males")]
    if (any(namesdf=="male")) males <- df[, which(namesdf=="male")]
    if (any(namesdf=="females")) females <- df[, which(namesdf=="females")]
    if (any(namesdf=="female")) females <- df[, which(namesdf=="female")]
    if (any(namesdf=="n")) N <- df[, which(namesdf=="n")]
    if (any(namesdf=="sexed")) N <- df[, which(namesdf=="sexed")]
    if (any(namesdf=="temperatures")) temperatures <- df[, which(namesdf=="temperatures")]
    if (any(namesdf=="temperature")) temperatures <- df[, which(namesdf=="temperature")]
    if (any(namesdf=="incubation.temperature")) temperatures <- df[, which(namesdf=="incubation.temperature")]
    if (any(namesdf=="durations")) durations <- df[, which(namesdf=="durations")]
    if (any(namesdf=="duration")) durations <- df[, which(namesdf=="duration")]
    if (any(namesdf=="IP.mean")) durations <- df[, which(namesdf=="IP.mean")]
  }
  
  if (is.null(durations) & is.null(temperatures)) {
    stop("Or temperatures or durations must be provided")
  }
  
  if (!is.null(durations) & !is.null(temperatures)) {
    stop("Temperatures and durations cannot be both provided")
  }
  
  
  # si je retire des lignes, c'est decale
  #  if (!is.null(males)) males <- na.omit(males)
  #  if (!is.null(females)) females <- na.omit(females)
  #  if (!is.null(N)) N <- na.omit(N)
  if (is.null(temperatures)) {
    IP <- TRUE
    temperatures <- durations
  } else {
    IP <- FALSE
  }
  
  cpt1 <- ifelse(is.null(males), 0, 1)+ifelse(is.null(females), 0, 1)+ifelse(is.null(N), 0, 1)
  if (cpt1<2) {
    stop("At least two informations from males, females or N must be provided")
  }
  
  if (is.null(males)) males <- N-females
  if (is.null(females)) females <- N-males
  if (is.null(N)) N <- males+females
  
  if (length(temperatures)!=length(N)) {
    stop("A temperature or duration must be provided for each experiment")
  }
  
  if (equation!="gsd" & equation!="hulin" & equation!="logistic" & equation!="hill" & equation!="richards" & equation!="double-richards") {
    stop("Equations supported are GSD, logistic, Hill, Hulin, Double_Richards and Richards")
  }
  
  
  if (equation=="gsd") {
    result <- list(par = NULL, SE=NULL, hessian=NULL, 
                   TRT=NULL,
                   SE_TRT=NULL,
                   message=NULL,
                   convergence=0,
                   counts=NULL,
                   value=-sum(dbinom(x=males, size=males+females, prob=0.5, log=TRUE)),
                   AIC=-2*sum(dbinom(x=males, size=males+females, prob=0.5, log=TRUE)))
#    limit.low.TRT <- min(temperatures)
#    limit.high.TRT <- max(temperatures)	
  } else {
    par <- parameters.initial
    if (!is.null(par)) {
    if (is.na(par["P"])) {
      par["P"] <- temperatures[which.min(abs((males/(males+females))-0.5))]
      if (IP) par["S"] <- abs(par["S"])
    }
      
    }
    if (equation!="richards") par <- par[which(names(par)!="K")]
    if (equation!="hulin" & equation!="double-richards") par <- par[which(names(par)!="K1")]
    if (equation!="hulin" & equation!="double-richards") par <- par[which(names(par)!="K2")]
    
    repeat {
     # result  <- optim(par, embryogrowth:::.tsd_fit, fixed.parameters=fixed.parameters, males=males, N=N, temperatures=temperatures, equation=equation, method="BFGS", hessian=TRUE, control = list(maxit=1000))
     result  <- optim(par, getFromNamespace(".tsd_fit", ns="embryogrowth"), 
           fixed.parameters=fixed.parameters, males=males, N=N, temperatures=temperatures, 
           equation=equation, method="BFGS", hessian=TRUE, control = list(maxit=1000))

      if (result$convergence != 1) break
      par<-result$par
      if (print) print("Convergence is not acheived. Optimization continues !")
    }

    par <- c(result$par, fixed.parameters)
    
    rh <- SEfromHessian(result$hessian, hessian=TRUE)
    result$SE <- rh$SE
    result$hessian <- rh$hessian
    result$AIC <- 2*result$value+2*length(par)
  }
  
  result$range.CI <- range.CI
  result$l <- l
  result$equation <- equation
  
  if (equation!="gsd") {

    class(result) <- "tsd"
   # result <<- result
    o <- P_TRT(x=result, temperatures=NULL, l=l, 
               replicate.CI = replicate.CI, probs = c((1-range.CI)/2, 0.5, 1-(1-range.CI)/2))
    
    result$P_TRT <- o$P_TRT_quantiles
  }
  
  saturated <- 0
  for (i in 1:length(N)) saturated <- saturated- dbinom(males[i], N[i], males[i]/N[i], log=TRUE)
  LRT <- 2*(result$value-saturated)
  result$GOF <- 1-pchisq(LRT, length(N)-length(par))
  
  result$males <- males
  result$females <- females
  result$N <- N
  result$temperatures <- temperatures
  result$males.freq <- males.freq

  result$fixed.parameters <- fixed.parameters
  # print(IP)
  result$type <- ifelse(IP, "duration", "temperature")
  
  if (print) print(paste("The goodness of fit test is", sprintf("%.5f",result$GOF)))
  
  if (equation!="gsd" & print) {
    if (!is.null(replicate.CI)) {
    print(paste0("The pivotal ", result$type, " is ", sprintf("%.3f",o$P_TRT_quantiles[2, "PT"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f", min(o$P_TRT_quantiles[1, "PT"], o$P_TRT_quantiles[3, "PT"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "PT"], o$P_TRT_quantiles[3, "PT"]))))
    print(paste0("The transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "TRT"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "TRT"], o$P_TRT_quantiles[3, "TRT"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "TRT"], o$P_TRT_quantiles[3, "TRT"]))))
    print(paste0("The lower limit of transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "lower.limit.TRT"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "lower.limit.TRT"], o$P_TRT_quantiles[3, "lower.limit.TRT"])), ";", sprintf("%.3f", max(o$P_TRT_quantiles[1, "lower.limit.TRT"], o$P_TRT_quantiles[3, "lower.limit.TRT"]))))
    print(paste0("The higher limit of transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "higher.limit.TRT"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "higher.limit.TRT"], o$P_TRT_quantiles[3, "higher.limit.TRT"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "higher.limit.TRT"], o$P_TRT_quantiles[3, "higher.limit.TRT"]))))
    } else {
      print(paste0("The pivotal ", result$type, " is ", sprintf("%.3f",o$P_TRT_quantiles[2, "PT"])))
      print(paste0("The transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "TRT"])))
      print(paste0("The lower limit of transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "lower.limit.TRT"])))
      print(paste0("The higher limit of transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "higher.limit.TRT"])))
    }
  }
  
  class(result) <- "tsd"
  return(invisible(result))
}
