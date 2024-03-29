#' tsd estimates the parameters that best describe temperature-dependent sex determination
#' @title Estimate the parameters that best describe temperature-dependent sex determination
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
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
#' @param equation Can be "logistic", "Hill", "A-logistic", "Hulin", "Double-A-logistic", "flexit", "GSD", "logit", "probit"
#' @param replicate.CI Number of replicates to estimate confidence intervals
#' @param replicate.NullDeviance Number of replicates to estimate null distribution of deviance
#' @param SE If FALSE, does not estimate SE of parameters. Can be use when something wrong happens.
#' @param range.CI The range of confidence interval for estimation, default=0.95
#' @param print Should the results be printed at screen? TRUE (default) or FALSE
#' @param control List of parameters used in optim.
#' @description Estimate the parameters that best describe the thermal reaction norm for sex ratio when temperature-dependent sex determination occurs.\cr
#' It can be used also to evaluate the relationship between incubation duration and sex ratio.\cr
#' The parameter l was defined in Girondot (1999). The TRT is defined from the difference between the two boundary temperatures giving sex ratios of \eqn{l} and \eqn{1 - l}, respectively:\cr
#' For logistic model (Girondot, 1999), it follows \deqn{TRT_{l}=abs\left ( S\: K_{l} \right )}{TRTl = abs( S.Kl )}
#' where \eqn{K_{l}}{Kl} is a constant equal to \eqn{2\: log\left ( \frac{l}{1-l} \right )}{2 log(l / ( 1 - l))}.\cr
#' In Girondot (1999), l was 0.05 and then the TRT was defined as being the range of temperatures producing from 5\% to 95\% of each sex.\cr
#' For other models, TRT is calculated numerically.\cr
#' The basic model is logistic one. This model has the particularity to have a symmetric shape around P.\cr
#' The other models have been built to alleviate this constraint. Hill and A-logistic models can be asymmetric, but it is impossible to control independently the low and high transitions.\cr
#' Hulin model is assymmetric but the control of asymmetry is difficult to manage.\cr
#' If asymmetric model is selected, it is always better to use flexit model.
#' \deqn{if dose < P then (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}{if dose < P then (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}
#' \deqn{if dose > P then 1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}{if dose > P then 1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}
#' with:\cr
#'      \deqn{S1 = S/((4/K1)*(2^(-K1))^(1/K1+1)*(2^K1-1))}{S1 = S/((4/K1)*(2^(-K1))^(1/K1+1)*(2^K1-1))}
#'      \deqn{S2 = S/((4/K2)*(2^(-K2))^(1/K2+1)*(2^K2-1))}{S2 = S/((4/K2)*(2^(-K2))^(1/K2+1)*(2^K2-1))}
#' @references Girondot, M. 1999. Statistical description of temperature-dependent sex determination using maximum likelihood. Evolutionary Ecology Research, 1, 479-486.
#' @references Godfrey, M.H., Delmas, V., Girondot, M., 2003. Assessment of patterns of temperature-dependent sex determination using maximum likelihood model selection. Ecoscience 10, 265-272.
#' @references Hulin, V., Delmas, V., Girondot, M., Godfrey, M.H., Guillon, J.-M., 2009. Temperature-dependent sex determination and global change: are some species at greater risk? Oecologia 160, 493-506.
#' @references Abreu-Grobois, F.A., Morales-Mérida, B.A., Hart, C.E., Guillon, J.-M., Godfrey, M.H., Navarro, E., Girondot, M., 2020. Recent advances on the estimation of the thermal reaction norm for sex ratios. PeerJ 8, e8451.
#' @family Functions for temperature-dependent sex determination
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' CC_AtlanticSW <- subset(DatabaseTSD, RMU=="Atlantic, SW" & 
#'                           Species=="Caretta caretta" & (!is.na(Sexed) & Sexed!=0) &
#'                           !is.na(Correction.factor))
#' tsdL <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="logistic", replicate.CI=NULL))
#' tsdH <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Hill", replicate.CI=NULL))
#' tsdR <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="A-logistic", replicate.CI=NULL))
#' tsdF <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Flexit", replicate.CI=NULL))
#' tsdDR <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Double-A-logistic", replicate.CI=NULL))
#' gsd <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="GSD", replicate.CI=NULL))
#' compare_AIC(Logistic_Model=tsdL, Hill_model=tsdH, Alogistic_model=tsdR, 
#'                flexit=tsdF, 
#'                DoubleAlogistic_model=tsdDR, GSD_model=gsd)
#' compare_AICc(Logistic_Model=tsdL, Hill_model=tsdH, Alogistic_model=tsdR, 
#'                DoubleAlogistic_model=tsdDR, GSD_model=gsd, factor.value = -1)
#' compare_BIC(Logistic_Model=tsdL, Hill_model=tsdH, Alogistic_model=tsdR, 
#'                DoubleAlogistic_model=tsdDR, GSD_model=gsd, factor.value = -1)
#' ##############
#' eo <- subset(DatabaseTSD, Species=="Emys orbicularis", c("Males", "Females", 
#'                                        "Incubation.temperature"))
#'                                        
#' eo_Hill <- with(eo, tsd(males=Males, females=Females, 
#'                                        temperatures=Incubation.temperature,
#'                                        equation="Hill"))
#' eo_Hill <- tsd(df=eo, equation="Hill", replicate.CI=NULL)
#' eo_logistic <- tsd(eo, replicate.CI=NULL)
#' eo_Alogistic <- with(eo, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature, 
#'                                  equation="a-logistic", replicate.CI=NULL))
#' ### The Hulin model is a modification of A-logistic (See Hulin et al. 2009)
#' 
#' ########## Caution
#' ### It should not be used anymore as it can produce unexpected results
#' par <- eo_Alogistic$par
#' names(par)[which(names(par)=="K")] <- "K2"
#' par <- c(par, K1=0)
#' eo_Hulin <- with(eo, tsd(males=Males, females=Females, 
#'                                  parameters.initial=par, 
#'                                  temperatures=Incubation.temperature, 
#'                                  equation="Hulin", replicate.CI=NULL))
#'                                  
#' ### The Double-A-logistic model is a A-logistic model with K1 and K2 using respectively
#' ### below and above P
#' 
#' ########## Caution
#' ### The curve is not smooth at pivotal temperature
#' 
#' par <- eo_Alogistic$par
#' names(par)[which(names(par)=="K")] <- "K2"
#' par <- c(par, K1=as.numeric(par["K2"])*0.8)
#' par["K2"] <- par["K2"]*0.8
#' eo_Double_Alogistic <- with(eo, tsd(males=Males, females=Females,
#'                                  parameters.initial=par,
#'                                  temperatures=Incubation.temperature,
#'                                  equation="Double-a-logistic", replicate.CI=NULL))
#'                                  
#' ### The flexit model is modeled with K1 and K2 using respectively
#' ### below and above P and smooth transition at P; S is the slope at P
#' 
#' par <- c(eo_logistic$par["P"], 1/4*eo_logistic$par["S"], K1=1, K2=1)
#' eo_flexit <- with(eo, tsd(males=Males, females=Females,
#'                                  parameters.initial=par,
#'                                  temperatures=Incubation.temperature,
#'                                  equation="flexit", replicate.CI=NULL))
#'                                  
#' compare_AIC(Logistic=eo_logistic, Hill=eo_Hill, Alogistic=eo_Alogistic, 
#'              Hulin=eo_Hulin, Double_Alogistic=eo_Double_Alogistic, 
#'              flexit=eo_flexit)
#' ## Note that SE for lower limit of TRT is wrong
#' plot(eo_flexit)
#' ## To get correct confidence interval, check \code{tsd_MHmcmc()}. 
#'              
#' ### Note the asymmetry of the Double-A-logistic and flexit models
#' predict(eo_Double_Alogistic, 
#'        temperatures=c(eo_Double_Alogistic$par["P"]-0.2, eo_Double_Alogistic$par["P"]+0.2))
#' predict(eo_Double_Alogistic)
#' 
#' (p <- predict(eo_flexit, 
#'        temperatures=c(eo_flexit$par["P"]-0.3, eo_flexit$par["P"]+0.3)))
#' p["50%", 1]-0.5; 0.5-p["50%", 2]
#' predict(eo_flexit)
#' 
#' ### It can be used also for incubation duration
#' CC_AtlanticSW <- subset(DatabaseTSD, RMU=="Atlantic, SW" & 
#'                           Species=="Caretta caretta" & Sexed!=0)
#' tsdL_IP <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  durations=IP.mean, 
#'                                  equation="logistic", replicate.CI=NULL))
#' plot(tsdL_IP, xlab="Incubation durations in days")
#' # Example with Chelonia mydas
#' cm <- subset(DatabaseTSD, Species=="Chelonia mydas" & !is.na(Sexed), c("Males", "Females", 
#'                                        "Incubation.temperature", "RMU"))
#' tsd(subset(cm, subset=RMU=="Pacific, SW"))
#' tsd(subset(cm, subset=RMU=="Pacific, Northwest"))
#' tsd(subset(cm, subset=RMU=="Atlantic, S Caribbean"))
#' 
#' ### Eretmochelys imbricata
#' Ei_PacificSW <- subset(DatabaseTSD, RMU=="Pacific, SW" & 
#'                        Species=="Eretmochelys imbricata")
#' Ei_AtlanticW <- subset(DatabaseTSD, RMU=="Atlantic, W (Caribbean and E USA)" & 
#'                        Species=="Eretmochelys imbricata")
#' Ei_AtlanticSW <- subset(DatabaseTSD, RMU=="Atlantic, SW" & 
#'                        Species=="Eretmochelys imbricata")
#' Ei_PacSW <- tsd(Ei_PacificSW)
#' Ei_AtlW <- tsd(Ei_AtlanticW)
#' Ei_AtlSW <- tsd(Ei_AtlanticSW)
#' 
#' plot(Ei_PacSW, xlim=c(27, 33), show.PTRT = FALSE, main=expression(italic("Eretmochelys imbricata")))
#' par(new=TRUE)
#' plot(Ei_AtlW, xlim=c(27, 33), col="red", xlab="", ylab="", 
#'      axes=FALSE, xaxt="n", show.PTRT = FALSE, errbar.col="red")
#' par(new=TRUE)
#' plot(Ei_AtlSW, xlim=c(27, 33), col="blue", xlab="", ylab="", axes=FALSE, 
#'      xaxt="n", show.PTRT = FALSE, errbar.col="blue")
#' legend("topright", legend=c("Pacific, SW", "Atlantic, W", "Atlantic, SW"), lty=1, 
#' col=c("black", "red", "blue"))
#' 
#' ### Chelonia mydas
#' Cm_PacificSW <- subset(DatabaseTSD, RMU=="Pacific, SW" & !is.na(Sexed) & 
#'                        Species=="Chelonia mydas")
#' Cm_PacificNW <- subset(DatabaseTSD, RMU=="Pacific, NW" &  !is.na(Sexed) & 
#'                        Species=="Chelonia mydas")
#' Cm_AtlanticSC <- subset(DatabaseTSD, RMU=="Atlantic, S Caribbean" &  !is.na(Sexed) & 
#'                        Species=="Chelonia mydas")
#' Cm_IndianSE <- subset(DatabaseTSD, RMU=="Indian, SE" &  !is.na(Sexed) & 
#'                        Species=="Chelonia mydas")
#' Cm_PacSW <- tsd(Cm_PacificSW)
#' Cm_PacNW <- tsd(Cm_PacificNW)
#' Cm_IndSE <- tsd(Cm_IndianSE)
#' Cm_AtlSC <- tsd(Cm_AtlanticSC)
#' 
#' plot(Cm_PacSW, xlim=c(24, 34), show.PTRT = FALSE, main=expression(italic("Chelonia mydas")))
#' par(new=TRUE)
#' plot(Cm_PacNW, xlim=c(24, 34), col="red", xlab="", ylab="", 
#'      axes=FALSE, xaxt="n", show.PTRT = FALSE, errbar.col="red")
#' par(new=TRUE)
#' plot(Cm_IndSE, xlim=c(24, 34), col="blue", xlab="", ylab="", 
#'      axes=FALSE, xaxt="n", show.PTRT = FALSE, errbar.col="blue")
#' par(new=TRUE)
#' plot(Cm_AtlSC, xlim=c(24, 34), col="green", xlab="", ylab="", 
#'      axes=FALSE, xaxt="n", show.PTRT = FALSE, errbar.col="green")
#' 
#' # To fit a TSDII or FMF TSD pattern, you must indicate P_low, S_low, P_high, and S_high
#' # for logistic model and P_low, S_low, K1_low, K2_low, P_high, S_high, K1_high, and K2_high for 
#' # flexit model
#' # The model must be 0-1 for low and 1-0 for high with P_low < P_high
#' 
#' Chelydra_serpentina <- subset(DatabaseTSD, !is.na(Sexed) & (Sexed != 0) & 
#'                        Species=="Chelydra serpentina")
#'                        
#' model_TSDII <- tsd(Chelydra_serpentina, males.freq=FALSE, 
#'                    parameters.initial=c(P_low=21, S_low=0.3, P_high=28, S_high=-0.4), 
#'                    equation="logistic")
#' plot(model_TSDII, lab.TRT = "TRT l = 5 %")
#' priors <- tsd_MHmcmc_p(result=model_TSDII, accept=TRUE)
#' out_mcmc <- tsd_MHmcmc(result=model_TSDII, n.iter=10000, parametersMCMC=priors)
#' plot(model_TSDII, resultmcmc=out_mcmc, lab.TRT = "TRT l = 5 %")
#' predict(model_TSDII, temperatures=25:35)
#' 
#' # Podocnemis expansa
#' Podocnemis_expansa <- subset(DatabaseTSD, !is.na(Sexed) & (Sexed != 0) & 
#'                        Species=="Podocnemis expansa")
#' Podocnemis_expansa_Valenzuela_2001 <- subset(Podocnemis_expansa, 
#'                    Reference=="Valenzuela, 2001")
#' PeL2001 <- tsd(df=Podocnemis_expansa_Valenzuela_2001)
#' # The pivotal temperature is 32.133 °C (CI 95% 31.495;32.766)
#' # In Valenzuela, 2001: "Using data from the present study alone, 
#' # the critical temperature was 32.2 °C by both methods and the 95% 
#' # confidence limits were 31.4 °C and 32.9 °C."
#' # Data are close but not identical to what was published.
#' 
#' # The pivotal temperature calculated by maximum likelihood and by inverse 
#' # prediction from logistic regression, was 32.6°C using raw data from 
#' # 1991 (N. Valenzuela, unpublished data) and from this study. The lower 
#' # and upper 95% confidence limits of the pivotal temperature were 32.2°C 
#' # and 33.2°C,
#' 
#' Podocnemis_expansa_Valenzuela_1997 <- subset(Podocnemis_expansa, 
#'                    subset=(((Reference=="Lance et al., 1992; Valenzuela et al., 1997") | 
#'                    (Reference=="Valenzuela, 2001")) & 
#'                    (!is.na(Sexed)) & (Sexed != 0)))
#' 
#' PeL1997 <- tsd(df=Podocnemis_expansa_Valenzuela_1997)
#' 
#' # Gekko japonicus
#' 
#' Gekko_japonicus <- subset(DatabaseTSD, !is.na(Sexed) & (Sexed != 0) & 
#' Species=="Gekko japonicus")
#' model_TSDII_gj <- tsd(Gekko_japonicus, males.freq=TRUE, 
#'                    parameters.initial=c(P_low=26, S_low=1.5, 
#'                                         P_high=31, S_high=-1.5), 
#'                    equation="logistic")
#' plot(model_TSDII_gj, lab.TRT = "TRT l = 5 %")
#' print(model_TSDII_gj)
#' prior <- tsd_MHmcmc_p(result = model_TSDII_gj, accept = TRUE)
#' prior <- structure(list(
#'   Density = c("dnorm", "dnorm", "dnorm", "dnorm"), 
#'   Prior1 = c(26, 0.3, 31, -0.4), 
#'   Prior2 = c(2, 1, 2, 1), 
#'   SDProp = c(2, 0.5, 2, 0.5), 
#'   Min = c(25, -2, 25, -2), 
#'   Max = c(35, 2, 35, 2), 
#'   Init = c(26, 0.3, 31, -0.4)), 
#'   row.names = c("P_low",  "S_low", "P_high", "S_high"), 
#'   class = "data.frame")
#'   
#' result_mcmc_tsd_gj <- tsd_MHmcmc(result=model_TSDII_gj, 
#' parametersMCMC=prior, n.iter=10000, n.chains = 1,  
#' n.adapt = 0, thin=1, trace=FALSE, adaptive=TRUE)
#' summary(result_mcmc_tsd_gj)
#' plot(result_mcmc_tsd_gj, parameters="P_low", scale.prior=TRUE, xlim=c(20, 30), las=1)
#' plot(result_mcmc_tsd_gj, parameters="P_high", scale.prior=TRUE, xlim=c(25, 35), las=1)
#' plot(model_TSDII_gj, resultmcmc = result_mcmc_tsd_gj)
#' 
#' # Trachemys scripta elegans
#' # Take care, the pattern reflects large population variation
#' 
#' Tse <- subset(DatabaseTSD, Species=="Trachemys scripta" & Subspecies == "elegans" & !is.na(Sexed))
#' Tse_logistic <- tsd(Tse)
#' plot(Tse_flexit)
#' compare_AICc(logistic=Tse_logistic, flexit=Tse_flexit)
#' plot(Tse_flexit)
#' }
#' @export


tsd <- function(df=NULL                                            , 
                males=NULL                                         , 
                females=NULL                                       , 
                N=NULL                                             , 
                temperatures=NULL                                  , 
                durations=NULL                                     ,
                l=0.05                                             , 
                parameters.initial=c(P=NA, S=-2, K=0, K1=1, K2=0)  , 
                males.freq=TRUE                                    , 
                fixed.parameters=NULL                              , 
                equation="logistic"                                , 
                replicate.CI=10000                                 , 
                range.CI=0.95                                      , 
                SE=TRUE                                            , 
                replicate.NullDeviance=1000                        , 
                control=list(maxit=1000)                           ,
                print=TRUE                                         ) {
  
  # df=NULL; males=NULL; females=NULL; N=NULL; temperatures=NULL; durations=NULL; l=0.05; parameters.initial=c(P=NA, S=-0.5, K=0, K1=1, K2=0); males.freq=TRUE; fixed.parameters=NULL; equation="logistic"; replicate.CI=10000; range.CI=0.95; SE=TRUE;print=TRUE; control=list(maxit=1000); replicate.NullDeviance=1000
  # CC_AtlanticSW <- subset(DatabaseTSD, RMU=="Atlantic, SW" & Species=="Caretta caretta" & Sexed!=0)
  # males=CC_AtlanticSW$Males; females=CC_AtlanticSW$Females; temperatures=CC_AtlanticSW$Incubation.temperature-CC_AtlanticSW$Correction.factor; equation="logistic"
  # parameters.initial = c(P=29.19); fixed.parameters = c(S=-0.21)
  equation <- tolower(equation)
  
  if (is.null(replicate.CI)) replicate.CI <- 0
  
  equation <- match.arg(equation, 
                        choices = c("logistic", "hill", "a-logistic", "hulin", 
                                    "double-a-logistic", "flexit", "gsd", 
                                    "probit", "logit"), 
                        several.ok = FALSE)
  
  if (!is.null(df)) {
    if ((!inherits(df, "data.frame")) & (!inherits(df, "matrix"))) {
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
  
  if (all(males * females == 0)) {
    warning("At least one incubation temperature should produce mixed sex ratio. Use MCMC rather.")
  }
  
  
  if (equation=="gsd") {
    value <- -sum(dbinom(x=males, size=males+females, prob=0.5, log=TRUE))
    result <- list(par = NULL, SE=NULL, hessian=NULL, 
                   TRT=NULL,
                   SE_TRT=NULL,
                   message=NULL,
                   convergence=0,
                   counts=NULL,
                   value = value,
                   AIC= 2*value + 2*0, 
                   AICc= 2*value + (2*0*(0+1))/(sum(N)-0-1), 
                   BIC= 2*value + 0
    )
    
    
    #    limit.low.TRT <- min(temperatures)
    #    limit.high.TRT <- max(temperatures)	
  } else {
    ppi <- parameters.initial
    par <- c(parameters.initial, fixed.parameters)
    if (is.na(par["P_low"])) {
      if (is.na(par["P"])) {
        if ((equation!="probit") & (equation!="logit")) {
          ppi["P"] <- temperatures[which.min(abs((males/(males+females))-0.5))]
        } else {
          ppi["P"] <- 0
        }
      }
      
      if (is.na(par["S"])) {
        if ((equation!="probit")) {
          pente <- lm(males/(males+females) ~ temperatures)
          ppi["S"] <- pente$coefficients["temperatures"]
          if (equation!="flexit") ppi["S"] <- 1/(4*par["S"])
        } else {
          ppi["S"] <- 0
        }
      }
      
    }
    
    if (IP) ppi["S"] <- abs(ppi["S"])
    
    if (equation != "a-logistic") ppi <- ppi[which(names(ppi)!="K")]
    if (equation != "hulin" & equation != "double-a-logistic" & equation != "flexit") {
      ppi <- ppi[which(names(ppi)!="K1")]
      ppi <- ppi[which(names(ppi)!="K2")]
    }
    
    par <-  ppi
    
    repeat {
      # result  <- optim(par, embryogrowth:::.tsd_fit, fixed.parameters=fixed.parameters, males=males, N=N, temperatures=temperatures, equation=equation, method="BFGS", hessian=TRUE, control = list(maxit=1000))
      result  <- optim(par, getFromNamespace(".tsd_fit", ns="embryogrowth"), 
                       fixed.parameters=fixed.parameters, males=males, N=N, temperatures=temperatures, 
                       equation=equation, method="BFGS", hessian=TRUE, control = control)
      
      if (result$convergence != 1) break
      par <- result$par
      if (print) print("Convergence is not acheived. Optimization continues !")
    }
    
    par <- c(result$par, fixed.parameters)
    
    if (SE) {
      rh <- SEfromHessian(result$hessian, hessian=TRUE)
      result$SE <- rh$SE
    }
    result$hessian <- result$hessian
    result$AIC <- 2*result$value+2*length(par)
    # Correction le 21/10/2020
    result$AICc <- result$AIC +(2*length(par)*(length(par)+1))/(sum(N)-length(par)-1)
    # Correction le 21/10/2020
    result$BIC <- 2*result$value+ length(par)*log(sum(N))
  }
  
  
  result$range.CI <- range.CI
  result$l <- l
  result$equation <- equation
  result$males <- males
  result$females <- females
  result$N <- N
  result$temperatures <- temperatures
  result$males.freq <- males.freq
  result$fixed.parameters <- fixed.parameters
  result$type <- ifelse(IP, "duration", "temperature")
  
  
  
  result$deviance <- -2*(-result$value - sum(dbinom(x = males, size=N, 
                                                    prob = males/N, log = TRUE)))
  # degrees of freedom calculated from the difference of the number of parameters in the saturated and the fitted model. 
  result$df <- length(males) - length(result$par)
  result$pvalue <- pchisq(q=result$deviance, df=result$df, lower.tail = FALSE)
  
  if (equation != "gsd") {
    
    result <- addS3Class(result, "tsd")
    # result <<- result
    o <- P_TRT(x=result, l=l, 
               replicate.CI = replicate.CI, 
               probs = c((1-range.CI)/2, 0.5, 1-(1-range.CI)/2))
    
    result$P_TRT <- o$P_TRT_quantiles
  }
  
  probtheor <- getFromNamespace(x=".modelTSD", ns="embryogrowth")(c(result$par, result$fixed.parameters), 
                                                                  result$temperatures, 
                                                                  equation=equation)
  result$predictMaleFrequency <- probtheor
  
  if (!is.null(replicate.NullDeviance) & (replicate.NullDeviance != 0)) {
    pt <- NULL
    # Bug correct 2020-01-11: change 1000 to replicate.NullDeviance
    for (i in 1:replicate.NullDeviance) {
      m <- rbinom(n=length(result$males), size=result$N, prob = probtheor)
      result_dev  <- optim(par=result$par, fn=getFromNamespace(".tsd_fit", ns="embryogrowth"), 
                           fixed.parameters=result$fixed.parameters, males=m, N=result$N, 
                           temperatures=result$temperatures, 
                           equation=result$equation, method="BFGS", 
                           hessian=FALSE, control = control)
      deviance_dev <- -2*(-result$value - sum(dbinom(x = m, size=N, 
                                                     prob = m/N, log = TRUE)))
      
      pt <- c(pt, deviance_dev)
    }
    
    result$NullDeviance <- pt
    result$NullDeviancePvalue <- sum(result$deviance < pt)/length(pt)
    
  }
  
  if (equation!="gsd" & print) {
    if (any(colnames(o$P_TRT_quantiles)=="PT")) {
      if (!is.null(replicate.CI) | (replicate.CI == 0)) {
        print(paste0("The pivotal ", result$type, " is ", sprintf("%.3f",o$P_TRT_quantiles[2, "PT"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f", min(o$P_TRT_quantiles[1, "PT"], o$P_TRT_quantiles[3, "PT"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "PT"], o$P_TRT_quantiles[3, "PT"]))))
        print(paste0("The transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "TRT"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "TRT"], o$P_TRT_quantiles[3, "TRT"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "TRT"], o$P_TRT_quantiles[3, "TRT"]))))
        print(paste0("The lower limit of transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "lower.limit.TRT"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "lower.limit.TRT"], o$P_TRT_quantiles[3, "lower.limit.TRT"])), ";", sprintf("%.3f", max(o$P_TRT_quantiles[1, "lower.limit.TRT"], o$P_TRT_quantiles[3, "lower.limit.TRT"]))))
        print(paste0("The upper limit of transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "higher.limit.TRT"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "higher.limit.TRT"], o$P_TRT_quantiles[3, "higher.limit.TRT"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "higher.limit.TRT"], o$P_TRT_quantiles[3, "higher.limit.TRT"]))))
        if (!is.na(result$par["S"])) print(paste0("The S parameter value is ", sprintf("%.3f",result$par["S"])))
      } else {
        print(paste0("The pivotal ", result$type, " is ", sprintf("%.3f",o$P_TRT_quantiles[2, "PT"])))
        print(paste0("The transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "TRT"])))
        print(paste0("The lower limit of transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "lower.limit.TRT"])))
        print(paste0("The upper limit of transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "higher.limit.TRT"])))
        if (!is.na(result$par["S"])) print(paste0("The S parameter value is ", sprintf("%.3f",result$par["S"])))
      }
    } else {
      if (!is.null(replicate.CI)) {
        print(paste0("The lower pivotal ", result$type, " is ", sprintf("%.3f",o$P_TRT_quantiles[2, "PT_low"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f", min(o$P_TRT_quantiles[1, "PT_low"], o$P_TRT_quantiles[3, "PT_low"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "PT_low"], o$P_TRT_quantiles[3, "PT_low"]))))
        print(paste0("The lower transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "TRT_low"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "TRT_low"], o$P_TRT_quantiles[3, "TRT_low"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "TRT_low"], o$P_TRT_quantiles[3, "TRT_low"]))))
        print(paste0("The lower limit of lower transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "lower.limit.TRT_low"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "lower.limit.TRT_low"], o$P_TRT_quantiles[3, "lower.limit.TRT_low"])), ";", sprintf("%.3f", max(o$P_TRT_quantiles[1, "lower.limit.TRT_low"], o$P_TRT_quantiles[3, "lower.limit.TRT_low"]))))
        print(paste0("The upper limit of lower transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "higher.limit.TRT_low"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "higher.limit.TRT_low"], o$P_TRT_quantiles[3, "higher.limit.TRT_low"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "higher.limit.TRT_low"], o$P_TRT_quantiles[3, "higher.limit.TRT_low"]))))
        print(paste0("The upper pivotal ", result$type, " is ", sprintf("%.3f",o$P_TRT_quantiles[2, "PT_high"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f", min(o$P_TRT_quantiles[1, "PT_high"], o$P_TRT_quantiles[3, "PT_high"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "PT_high"], o$P_TRT_quantiles[3, "PT_high"]))))
        print(paste0("The upper transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "TRT_high"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "TRT_high"], o$P_TRT_quantiles[3, "TRT_high"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "TRT_high"], o$P_TRT_quantiles[3, "TRT_high"]))))
        print(paste0("The lower limit of upper transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "lower.limit.TRT_high"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "lower.limit.TRT_high"], o$P_TRT_quantiles[3, "lower.limit.TRT_high"])), ";", sprintf("%.3f", max(o$P_TRT_quantiles[1, "lower.limit.TRT_high"], o$P_TRT_quantiles[3, "lower.limit.TRT_high"]))))
        print(paste0("The upper limit of upper transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "higher.limit.TRT_high"]), " CI", as.character(range.CI*100), "% ", sprintf("%.3f",min(o$P_TRT_quantiles[1, "higher.limit.TRT_high"], o$P_TRT_quantiles[3, "higher.limit.TRT_high"])), ";", sprintf("%.3f",max(o$P_TRT_quantiles[1, "higher.limit.TRT_high"], o$P_TRT_quantiles[3, "higher.limit.TRT_high"]))))
        if (!is.na(result$par["S"])) print(paste0("The S parameter value is ", sprintf("%.3f",result$par["S"])))
      } else {
        print(paste0("The lower pivotal ", result$type, " is ", sprintf("%.3f",o$P_TRT_quantiles[2, "PT_low"])))
        print(paste0("The lower transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "TRT_low"])))
        print(paste0("The lower limit of lower transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "lower.limit.TRT_low"])))
        print(paste0("The upper limit of lower transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "higher.limit.TRT_low"])))
        print(paste0("The upper pivotal ", result$type, " is ", sprintf("%.3f",o$P_TRT_quantiles[2, "PT_high"])))
        print(paste0("The upper transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "TRT_high"])))
        print(paste0("The lower limit of upper transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "lower.limit.TRT_high"])))
        print(paste0("The upper limit of upper transitional range of ", result$type, "s is ", sprintf("%.3f",o$P_TRT_quantiles[2, "higher.limit.TRT_high"])))
        if (!is.na(result$par["S"])) print(paste0("The S parameter value is ", sprintf("%.3f",result$par["S"])))
      }
    }
  }
  
  result <- addS3Class(result, "tsd")
  return(invisible(result))
}
