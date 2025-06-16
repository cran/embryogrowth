#' searchR fits the parameters that best represent nest incubation data.
#' @title Fit the parameters that best represent nest incubation data.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A result object
#' @param parameters A set of parameters used as initial point for searching
#' @param fixed.parameters A set of parameters that will not be changed
#' @param temperatures Timeseries of temperatures after formated using FormatNests()
#' @param integral Function used to fit embryo growth: integral.Gompertz, integral.exponential or integral.linear
#' @param derivate Function used to fit embryo growth: derivate.Gompertz, derivate.exponential or derivate.linear
#' @param hatchling.metric A vector with Mean and SD of size of hatchlings, ex. hatchling.metric=c(Mean=39, SD=3). Can be a data.frame also. See description
#' @param M0 Measure of hatchling size or mass proxi at laying date
#' @param saveAtMaxiter If True, each time number of interation reach maxiter, current data are saved in file with filename name
#' @param fileName The intermediate results are saved in file with fileName.Rdata name
#' @param weight A named vector of the weight for each nest for likelihood estimation
#' @param control List for control parameters for optimx
#' @description Fit the parameters that best represent data.\cr
#' hatchling.metric can be a data.frame with two columns Mean and SD and rownames with the nest name.\cr
#' If SD is na, then least squarre criteria is used for fitting.\cr
#' Function to fit thermal reaction norm can be expressed as : \cr
#' \itemize{
#'   \item a 4-parameters Schoolfield, Sharpe, and Magnuson model (1981) with \code{DHH}, \code{DHA}, \code{T12H}, and \code{Rho25};
#'   \item a 6-parameters Schoolfield, Sharpe, and Magnuson model (1981) with \code{T12L}, \code{DT}, \code{DHH}, \code{DHL}, \code{DHA}, and \code{Rho25};
#'   \item Each of these two first models can be combined as low and high sets of parameters by adding the \code{_L} suffix to one set. Then you must add also \code{transition_S} and \code{transition_P} parameters and then the growth rate is 1/(1+exp((1/transition_S)*(P-transition_P))) with P being the proportion of development;
#'   \item The \code{Rho25_b} control the effect of hygrometry (or \code{Rho25_b_L}) (It is not fully functional still);
#'   \item a Weibull function with \code{k} (shape), \code{lambda} (scale) and \code{theta} parameters;
#'   \item a normal function with \code{Peak}, \code{Scale}, and \code{sd} parameters;
#'   \item an asymmetric normal fuction with \code{Peak}, \code{Scale}, \code{sdH} and \code{sdL} parameters;
#'   \item a symmetric trigonometric function with \code{Length}, \code{Peak}, and \code{Max};
#'   \item an asymmetric trigonometric function with \code{LengthB}, \code{LengthE}, \code{Peak}, and \code{Max}.
#'   \item Dallwitz-Higgins model (1992) can be used using \code{Dallwitz_b1}, \code{Dallwitz_b2}, \code{Dallwitz_b3}, \code{Dallwitz_b4} and \code{Dallwitz_b5} parameters.
#'   \item If \code{Dallwitz_b4} is not included, \code{Dallwitz_b4} = 6 will be used.
#'   \item If \code{Dallwitz_b5} is not included, \code{Dallwitz_b5} = 0.4 will be used.
#'   \item It is possible also to add the parameter \code{epsilon} and then the model becomes \code{X + epsilon} with X being any of the above model;
#'   \item It is possible also to add the parameter \code{epsilon_L} and then the model becomes \code{X_L + epsilon_L} with \code{X_L} being any of the above model with suffix \code{_L};
#'   \item If the name of the parameter is a number, then the model is a polynom anchored with the rate being the parameter value at this temperature (the name). see \code{ChangeSSM()} function.
#' }
#' @references
#' \insertRef{9039}{embryogrowth}\cr
#' \insertRef{10871}{embryogrowth}\cr
#' \insertRef{8566}{embryogrowth}\cr
#' \insertRef{10620}{embryogrowth}\cr
#' \insertRef{10039}{embryogrowth}\cr
#' \insertRef{3326}{embryogrowth}\cr
#' \insertRef{5322}{embryogrowth}\cr
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' Laying.Time <- matrix(c("DY.1", "15/05/2010", 
#'                  "DY.17", "24/05/2010", 
#'                  "DY.16", "24/05/2010", 
#'                  "DY.18", "25/05/2010", 
#'                  "DY.20", "25/05/2010", 
#'                  "DY.21", "26/05/2010", 
#'                  "DY.22", "26/05/2010", 
#'                  "DY.23", "26/05/2010", 
#'                  "DY.24", "27/05/2010", 
#'                  "DY.25", "27/05/2010", 
#'                  "DY.28", "28/05/2010", 
#'                  "DY.26", "28/05/2010", 
#'                  "DY.27", "28/05/2010", 
#'                  "DY.146", "20/06/2010", 
#'                  "DY.147", "20/06/2010", 
#'                  "DY.172", "24/06/2010", 
#'                  "DY.175", "24/06/2010", 
#'                  "DY.170", "24/06/2010", 
#'                  "DY.260", "06/07/2010", 
#'                  "DY.282", "12/07/2010", 
#'                  "DY.310", "18/07/2010", 
#'                  "DY.309", "18/07/2010", 
#'                  "DY.328", "25/07/2010", 
#'                  "DY.331", "26/07/2010"), byrow=TRUE, ncol=2)
#' tz <- OlsonNames()[grepl("Asia/Istanbul", OlsonNames())]
#' Laying.Time_f <- setNames(as.POSIXlt.character(Laying.Time[, 2], format = "%d/%m/%Y", tz=tz), 
#'                            Laying.Time[, 1])
#' formated <- FormatNests(data=nest, previous=NULL, col.Time="Time", 
#'                         LayingTime=Laying.Time_f, 
#'                         hatchling.metric.mean=39.33, 
#'                         hatchling.metric.sd=1.92)
#' plot(formated, series=c(1, 2), lwd=3)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' # K for Gompertz must be set as fixed parameter or being a constant K  
#' # or relative to the hatchling size rK
#' ############################################################################
#' # From Girondot M, Monsinjon J, Guillon J-M (2018) Delimitation of the 
#' # embryonic thermosensitive period for sex determination using an embryo 
#' # growth model reveals a potential bias for sex ratio prediction in turtles. 
#' # Journal of Thermal Biology 73: 32-40 
#' #  rK = 1.208968 
#' #  M0 = 0.3470893 
#' ############################################################################
#' pfixed <- c(rK=1.208968)
#' M0 = 0.3470893 
#' ############################################################################
#' # 4 parameters SSM
#' ############################################################################
#' x4 <- c('DHA' = 109.31113503282113, 'DHH' = 617.80695919563857, 
#'     'T12H' = 306.38890489505093, 'Rho25' = 229.37265815800225)
#'     
#' resultNest_4p_SSM <- searchR(parameters=x4, fixed.parameters=pfixed, 
#' 	                            temperatures=formated, 
#' 	                            integral=integral.Gompertz, M0=M0)
#' 	
#' data(resultNest_4p_SSM)
#' plot(resultNest_4p_SSM, xlim=c(0,70), ylimT=c(22, 32), 
#'      ylimS=c(0,45), series=1, 
#'      embryo.stages="Caretta caretta.SCL")
#' plotR(resultNest_4p_SSM, ylim=c(0,6))
#' 
#' ############################################################################
#' # 6 parameters SSM
#' ############################################################################
#' x6 <- structure(c(106.567809092008, 527.359011254683, 614.208632495199, 
#' 2720.94506457237, 306.268259715624, 120.336791245212), .Names = c("DHA", 
#' "DHH", "DHL", "DT", "T12L", "Rho25"))
#' 
#' ############################################################################
#' # example of data.frame for hatchling.metric
#' ############################################################################
#' thatchling.metric <- data.frame(Mean = rep(39.33, formated$IndiceT["NbTS"]), 
#'                      SD = rep(1.92, formated$IndiceT["NbTS"]), 
#'                      row.names = formated$Names)
#' # It is sometimes difficult to find a good starting point for 
#' # SSM 6 parameters model. This function helps to find it based on a previoulsy 
#' # fitted model.
#' 
#'  x <- ChangeSSM(temperatures = (200:350)/10,
#'                 parameters = resultNest_4p_SSM$par,
#'                 initial.parameters = x6, 
#'                 control=list(maxit=1000))
#'                      
#' resultNest_6p_SSM <- searchR(parameters=x$par, fixed.parameters=pfixed, 
#' 	                        temperatures=formated, integral=integral.Gompertz, 
#' 	                        M0=M0, 
#' 	                        hatchling.metric=thatchling.metric)
#' 	                        
#' plotR(resultNest_6p_SSM, curve = "ML", ylim=c(0, 8))
#' 	                        
#' data(resultNest_6p_SSM)
#' pMCMC <- TRN_MHmcmc_p(resultNest_6p_SSM, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' resultNest_mcmc_6p_SSM <- GRTRN_MHmcmc(result=resultNest_6p_SSM,  
#' 	                                      parametersMCMC=pMCMC, 
#' 	                                      n.iter=10000, 
#' 	                                      n.chains = 1, 
#' 	                                      n.adapt = 0, 
#' 	                                      thin=1, 
#' 	                                      trace=TRUE)
#' 	
#' ############################################################################
#' # compare_AIC() is a function from the package "HelpersMG"
#' compare_AIC(test1=resultNest_4p_SSM, test2=resultNest_6p_SSM)
#' ############################################################################
#' 
#' ############################################################################
#' ############ Example as linear progression of development
#' ############ The development progress goes from 0 to 1
#' ############################################################################
#' 
#' pfixed <- NULL
#' M0 = 0
#' 
#' ############################################################################
#' # 4 parameters SSM
#' ############################################################################
#' x4 <- c('DHA' = 64.868697530424186, 'DHH' = 673.18292743646771, 
#'        'T12H' = 400.90952554047749, 'Rho25' = 82.217237723502123)
#' resultNest_4p_SSM_Linear <- searchR(parameters=x4, fixed.parameters=pfixed, 
#' 	temperatures=formated, integral=integral.linear, M0=M0, 
#' 	hatchling.metric=c(Mean=39.33, SD=1.92)/39.33)
#' plotR(resultNest_4p_SSM_Linear, ylim=c(0, 2), scaleY= 100000, curve = "ML")
#' plot(resultNest_4p_SSM_Linear, xlim=c(0,70), ylimT=c(22, 32), ylimS=c(0,1.1), 
#' series=1, embryo.stages="Generic.ProportionDevelopment")
#' 
#' tc <- GenerateConstInc(duration=600*24*60, temperatures = 28)
#' tc_f <- FormatNests(tc)
#' 
#' plot(x=resultNest_4p_SSM_Linear, xlim=c(0,70), ylimT=c(22, 32), ylimS=c(0,1.1), 
#'      series=1, embryo.stages="Generic.ProportionDevelopment", 
#'      stop.at.hatchling.metric=TRUE, metric.end.incubation="hatchling.metric", 
#'      temperatures=tc_f, hatchling.metric=c(Mean=39.33, SD=1.92)/39.33, 
#'      show.TSP=FALSE)
#'
#' ############################################################################
#' ############ with new parametrization based on anchor
#' ############ This is a non-parametric version
#' ############################################################################
#' 
#' data(resultNest_4p_SSM)
#' x0 <- resultNest_4p_SSM$par
#' t <- range(hist(resultNest_4p_SSM, plot=FALSE)$temperatures)
#' 
#' x <- getFromNamespace(".SSM", ns="embryogrowth")(T=seq(from=t[1], 
#'                                                        to=t[2], 
#'                                                        length.out=7), 
#'                                             parms=x0)[[1]]*1E5
#' names(x) <- as.character(seq(from=t[1], 
#'                           to=t[2], 
#'                           length.out=7))
#'      
#' M0 <- 0.3470893 
#' pfixed <- c(rK=1.208968)
#' resultNest_newp <- searchR(parameters=x, fixed.parameters=pfixed,
#'                            temperatures=formated, 
#'                            integral=integral.Gompertz, M0=M0, 
#'                            hatchling.metric=c(Mean=39.33, SD=1.92))
#' plotR(resultNest_newp, ylim=c(0, 2), 
#'       xlim=c(23, 34), ylimH=c(0, 3), show.hist=TRUE)
#' compare_AIC(test4p=resultNest_4p_SSM, 
#'             test6p=resultNest_6p_SSM, 
#'             testAnchor=resultNest_newp)
#'             
#' ############################################
#' # example with thermal reaction norm fitted from Weibull function
#' ############################################
#' 
#'  x <- ChangeSSM(temperatures = (200:350)/10,
#'                 parameters = resultNest_4p_SSM$par,
#'                 initial.parameters = structure(c(73.4009010417375, 304.142079511996, 
#'                                                 27.4671689276281), 
#'                                         .Names = c("k", "lambda", "scale")), 
#'                 control=list(maxit=1000))
#' M0 <- 0.3470893 
#' pfixed <- c(rK=1.208968)
#' resultNest_3p_Weibull <- searchR(parameters=x$par, fixed.parameters=pfixed, 
#'                          temperatures=formated, integral=integral.Gompertz, M0=M0, 
#'                          hatchling.metric=c(Mean=39.33, SD=1.92))
#' plotR(resultNest_3p_Weibull, ylim=c(0,6), col="Black")
#' compare_AIC(SSM=resultNest_4p_SSM, Weibull=resultNest_3p_Weibull)
#' 
#' ###########################################
#' # example with thermal reaction norm fitted from asymmetric normal function
#' ############################################
#' 
#' x <- ChangeSSM(temperatures = (200:350)/10,
#'                parameters = resultNest_4p_SSM$par,
#'                initial.parameters = structure(c(3, 7, 11, 32), 
#'                                .Names = c("Scale", "sdL", "sdH", "Peak")), 
#'                control=list(maxit=1000))
#' M0 <- 0.3470893 
#' pfixed <- c(rK=1.208968)
#' resultNest_4p_normal <- searchR(parameters=x$par, fixed.parameters=pfixed, 
#'                          temperatures=formated, integral=integral.Gompertz, M0=M0, 
#'                          hatchling.metric=c(Mean=39.33, SD=1.92))
#'                          
#' ###########################################
#' # example with thermal reaction norm fitted from trigonometric model
#' ############################################
#' 
#'  x <- ChangeSSM(temperatures = (200:350)/10,
#'                parameters = resultNest_4p_SSM$par,
#'                initial.parameters = structure(c(3, 20, 40, 32), 
#'                .Names = c("Max", "LengthB", "LengthE", "Peak")), 
#'                control=list(maxit=1000))
#' M0 <- 0.3470893 
#' pfixed <- c(rK=1.208968)
#' resultNest_4p_trigo <- searchR(parameters=x$par, fixed.parameters=pfixed, 
#'                          temperatures=formated, integral=integral.Gompertz, M0=M0, 
#'                          hatchling.metric=c(Mean=39.33, SD=1.92))
#'                          
#' ###############################################################
#' # Example with thermal reaction norm fitted from Dallwitz model
#' ###############################################################
#' # See: Dallwitz, M.J., Higgins, J.P., 1992. User’s guide to DEVAR. A computer 
#' # program for estimating development rate as a function of temperature. CSIRO Aust 
#' # Div Entomol Rep 2, 1-23.
#' 
#' # Note that Dallwitz model has many problems and I recommend to not use it:
#' # - The 3-parameters is too highly constraint
#' # - The 5 parameters produced infinite outputs for some sets of parameters that
#' #   can be generated while using delta method.
#' 
#' x <- c('Dallwitz_b1' = 4.8854060791241816, 
#'        'Dallwitz_b2' = 20.398366565842029, 
#'        'Dallwitz_b3' = 31.510995256647092)
#' M0 <- 0.3470893 
#' pfixed <- c(rK=1.208968)
#' resultNest_3p_Dallwitz <- searchR(parameters=x, fixed.parameters=pfixed, 
#'                          temperatures=formated, integral=integral.Gompertz, M0=M0, 
#'                          hatchling.metric=c(Mean=39.33, SD=1.92))
#' plotR(resultNest_3p_Dallwitz, ylim=c(0,6))
#' 
#' x <- c('Dallwitz_b1' = 4.9104386262684656, 
#'        'Dallwitz_b2' = 7.515425231891359, 
#'        'Dallwitz_b3' = 31.221784599026638, 
#'        'Dallwitz_b4' = 7.0354467023505682, 
#'        'Dallwitz_b5' = -1.5955717975708577)
#' pfixed <- c(rK=1.208968)
#' resultNest_5p_Dallwitz <- searchR(parameters=x, fixed.parameters=pfixed, 
#'                          temperatures=formated, integral=integral.Gompertz, M0=0.3470893, 
#'                          hatchling.metric=c(Mean=39.33, SD=1.92))
#' plotR(resultNest_5p_Dallwitz, ylim=c(0,3), scaleY=10000)
#' 
#' xp <- resultNest_6p_SSM$par
#' xp["Rho25"] <- 233
#' pfixed <- c(rK=1.208968)
#' resultNest_6p_SSM <- searchR(parameters=xp, fixed.parameters=pfixed, 
#'                          temperatures=formated, integral=integral.Gompertz, M0=0.3470893, 
#'                          hatchling.metric=c(Mean=39.33, SD=1.92))
#' plotR(resultNest_6p_SSM, ylim=c(0,8))
#' 
#' xp <- ChangeSSM(parameters = resultNest_3p_Dallwitz$par, 
#'                 initial.parameters = resultNest_4p_SSM$par)
#' pfixed <- c(rK=1.208968)
#' resultNest_4p_SSM <- searchR(parameters=xp$par, fixed.parameters=pfixed, 
#'                          temperatures=formated, integral=integral.Gompertz, M0=0.3470893, 
#'                          hatchling.metric=c(Mean=39.33, SD=1.92))
#' plotR(resultNest_4p_SSM, ylim=c(0,6))
#' 
#' compare_AIC(Dallwitz3p=resultNest_3p_Dallwitz, Dallwitz5p=resultNest_5p_Dallwitz, 
#'              SSM=resultNest_4p_SSM, SSM=resultNest_6p_SSM)
#'                          
#' ###########################################
#' # Example with thermal reaction norm of proportion of development
#' # fitted from Dallwitz model
#' # see Woolgar, L., Trocini, S., Mitchell, N., 2013. Key parameters describing 
#' # temperature-dependent sex determination in the southernmost population of loggerhead 
#' # sea turtles. Journal of Experimental Marine Biology and Ecology 449, 77-84.
#' ############################################
#' 
#' x <- structure(c(1.48207559695689, 20.1100310234046, 31.5665036287242), 
#'  .Names = c("Dallwitz_b1", "Dallwitz_b2", "Dallwitz_b3"))
#' resultNest_PropDev_3p_Dallwitz <- searchR(parameters=x, fixed.parameters=NULL, 
#'                          temperatures=formated, integral=integral.linear, M0=0, 
#'                          hatchling.metric=c(Mean=1, SD=NA))
#'  plotR(resultNest_PropDev_3p_Dallwitz, ylim=c(0, 1.5), curve="ML")
#'  plot(x=resultNest_PropDev_3p_Dallwitz, ylimS=c(0,1), xlim=c(0,60), series=2, 
#'          embryo.stages="Generic.ProportionDevelopment")
#'          
#' x <- structure(c(1.48904182113431, 10.4170365155993, 31.2591665490154, 
#' 6.32355497589913, -1.07425378667104), .Names = c("Dallwitz_b1", 
#' "Dallwitz_b2", "Dallwitz_b3", "Dallwitz_b4", "Dallwitz_b5"))
#' resultNest_PropDev_5p_Dallwitz <- searchR(parameters=x, fixed.parameters=NULL, 
#'                          temperatures=formated, integral=integral.linear, M0=0, 
#'                          hatchling.metric=c(Mean=1, SD=NA))
#'  plotR(resultNest_PropDev_5p_Dallwitz, ylim=c(0, 1.5))
#'  plot(x=resultNest_PropDev_5p_Dallwitz, ylimS=c(0,1), xlim=c(0,60), series=2, 
#'          embryo.stages="Generic.ProportionDevelopment")
#'          
#'  plotR(resultNest_PropDev_3p_Dallwitz, ylim=c(0, 1.5), curve="ML")
#'  plotR(resultNest_PropDev_5p_Dallwitz, ylim=c(0, 1.5), curve="ML", new=FALSE, col="red")
#'  compare_AICc(Dallwitz3p=resultNest_PropDev_3p_Dallwitz, 
#'               Dallwitz5p=resultNest_PropDev_5p_Dallwitz)
#'               
#' ###########################################################################
#' # Dalwitz model with proportion of development and fitted SD for final size
#' ###########################################################################
#' 
#' x <- c('Dallwitz_b1' = 1.4886497996404355, 
#'        'Dallwitz_b2' = 10.898310418085916, 
#'        'Dallwitz_b3' = 31.263224721068056, 
#'        'Dallwitz_b4' = 6.1624623077734535, 
#'        'Dallwitz_b5' = -1.0027132357973265, 
#'        'SD' = 0.041829475961912894)
#' resultNest_PropDev_5p_Dallwitz <- searchR(parameters=x, fixed.parameters=NULL, 
#'                          temperatures=formated, integral=integral.linear, M0=0, 
#'                          hatchling.metric=c(Mean=1))
#'  plotR(resultNest_PropDev_5p_Dallwitz, ylim=c(0, 1.5), curve="ML")
#'  # Note that the standard error of the curve cannot be estimated with delta method. 
#'  # MCMC should be used
#'  plot(x=resultNest_PropDev_5p_Dallwitz, ylimS=c(0,1), xlim=c(0,60), series=2, 
#'          embryo.stages="Generic.ProportionDevelopment")
#'          
#' ##############################################################################         
#' # Parameters Threshold_Low and Threshold_High are used to truncate growth rate
#' ##############################################################################         
#' 
#' plotR(result=resultNest_PropDev_5p_Dallwitz, 
#'      fixed.parameters=c(Threshold_Low=26, 
#'                         Threshold_High=33), 
#'      ylim=c(0, 1.5), curve="ML")
#' }
#' @export


searchR <- function(parameters=stop('Initial set of parameters must be provided'), 
                    temperatures=stop('Formated temperature must be provided')   , 
                    fixed.parameters=c(rK=1.208968)                              , 
                    integral=integral.Gompertz                                   , 
                    derivate=NULL                                                , 
                    hatchling.metric=c(Mean=39.33, SD=1.92)                      , 
                    M0=0.3470893                                                 , 
                    saveAtMaxiter=FALSE                                          , 
                    fileName="intermediate"                                      , 
                    weight=NULL                                                  , 
                    control=list(trace=1, REPORT=100, maxit=500)                 ) {
  
  # parameters=NULL; fixed.parameters=NULL; temperatures=NULL; integral=integral.Gompertz; derivate=NULL; hatchling.metric=c(Mean=39.33, SD=1.92); M0=1.7; saveAtMaxiter=FALSE; fileName="intermediate"; weight=NULL; SE=FALSE; control=list(trace=1, REPORT=100, maxit=500)  # temperatures <- formated
  
  # parameters=x; fixed.parameters=pfixed; temperatures=formated; integral=integral.Gompertz; M0=1.7; hatchling.metric=c(Mean=39.33, SD=1.92)
  
  mc.cores <- getOption("mc.cores", parallel::detectCores())
  message(paste0("I will use ", as.character(mc.cores)," cores."))
  
  if (!requireNamespace("optimx", quietly = TRUE)) {
    stop("optimx package is absent; Please install it first")
  }
  # library("optimx")
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("numDeriv package is absent; Please install it first")
  }
  # library("numDeriv")
  # dans temperatures il faut que je rajoute une colonne avec les indices de temperatures en K
  
  method <- c("Nelder-Mead","BFGS")
  NbTS <- temperatures$IndiceT["NbTS"]
  
  # si j'ai weight dans les data formatees et pas en paramtres, je les prends
  if (is.null(weight)) {
    weight <- unlist(lapply(temperatures$Nests, FUN = function(x) unname(x$weight)))
    if (is.null(weight)) {
      weight <- setNames(rep(1, NbTS), temperatures$Names)
    }
  }
  
  if (!setequal(names(weight), temperatures$Names)) {
    stop(paste("The weight parameter must define weight for each nest. Check", 
               setdiff(temperatures$Names, names(weight)), "nests/n"))
  }
  
  ##########################################################
  # Donnees de base de Gompertz
  ##########################################################
  
  if (is.null(hatchling.metric)) {
    hatchling.metric.mean <- unlist(lapply(temperatures$Nests, FUN = function(x) x$hatchling.metric.mean))
    hatchling.metric.sd <- unlist(lapply(temperatures$Nests, FUN = function(x) x$hatchling.metric.sd))
    hatchling.metricuse <- data.frame(Mean=hatchling.metric.mean, 
                                      SD=hatchling.metric.sd, 
                                      row.names=temperatures$Names)
  } else {
    if (is.numeric(hatchling.metric)) {
      hatchling.metricuse <- data.frame(Mean=rep(hatchling.metric["Mean"], NbTS), 
                                        SD=rep(hatchling.metric["SD"], NbTS), 
                                        row.names=temperatures$Names)
    } else {
      hatchling.metricuse <- hatchling.metric
    }
  }
  
  
  if (any(is.na(hatchling.metricuse[, "SD"]))) {
    hatchling.metricuse[, "SD"] <- NA
    if (is.na(c(fixed.parameters, parameters)["SD"])) {
      message("Fit criteria will be SSE.")
    } else {
      # 28/3/2022
      message("Fit criteria will be -Ln L with fitted sd for hatchling size.")
    }
  } else {
    message("Fit criteria will be -Ln L.")
  }
  
  # 25/2/2015
  for (j in temperatures$Names) temperatures$Nests[[j]]$data[1, "Mass"] <- M0
  
  # Un paramtre ne peut pas tre indique en fixe et en fite - 22/7/2012	
  
  if (length(intersect(names(parameters), names(fixed.parameters)))!=0) {
    stop("A parameter cannot be fixed and fitted at the same time !")
  }
  
  grR <- getFromNamespace(".gradientRichardson", ns="embryogrowth")
  nm <- names(parameters)
  
  repeat {
    result <- try(optimx::optimx(par=parameters, 
                                 fn=getFromNamespace("info.nests", ns="embryogrowth"), 
                                 temperatures=temperatures, 
                                 integral=integral, 
                                 weight=weight, 
                                 derivate=derivate, 
                                 WAIC=FALSE, 
                                 hatchling.metric=hatchling.metricuse, 
                                 M0=M0, 
                                 fixed.parameters=fixed.parameters,
                                 gr=grR, 
                                 method=method, 
                                 control=modifyList(list(dowarn=FALSE, follow.on=TRUE, kkt=FALSE), control), 
                                 hessian=FALSE), silent=TRUE)
    
    # getFromNamespace("info.nests", ns="embryogrowth")(parameters=parameters, 
    #                 temperatures=temperatures, 
    #                 integral=integral, weight=weight, derivate=derivate, 
    #                  hatchling.metric=hatchling.metricuse, M0=M0, 
    #                  fixed.parameters=fixed.parameters)
    
    # parameters=parameters
    # temperatures=temperatures
    # integral=integral; weight=weight; derivate=derivate, 
    # hatchling.metric=hatchling.metricuse; M0=M0, 
    # fixed.parameters=fixed.parameters
    
    minL <- nrow(result)
    
    x <- unlist(result[minL, 1:length(nm)])
    names(x) <- nm
    conv <- result[minL, "convcode"]
    value <- result[minL, "value"]
    
    if (conv != 1) break
    parameters <- x
    message("Convergence is not achieved. Optimization continues !")
    message(dput(parameters))
    if (saveAtMaxiter) save(parameters, file=paste0(fileName, ".RData"))
  }
  
  result_list <- list()
  result_list$result <- result
  result_list$par <- x
  result_list$value <- value
  result_list$convergence <- conv
  result <- result_list
  
  print(d(result$par))
  
  
  if ((any(is.na(hatchling.metricuse[, "SD"]))) & (is.na(x["SD"]))) {
    result$criteria <- "SSE"
    result$hessian <- NULL
    res <- rep(NA, length(parameters))
    names(res) <- names(parameters)
    
  } else {
    result$criteria <- "-Ln L"
    
    mathessian <- try(numDeriv::hessian(getFromNamespace("info.nests", ns="embryogrowth"), 
                                        x=result$par, 
                                        method="Richardson", 
                                        temperatures=temperatures, 
                                        integral=integral, 
                                        derivate=derivate, 
                                        weight=weight,
                                        hatchling.metric=hatchling.metricuse, 
                                        M0=M0, 
                                        fixed.parameters=fixed.parameters)
                      , silent=TRUE)
    
    if (inherits(mathessian, "try-error")) {
      result$hessian <- NULL
      res <- rep(NA, length(parameters))
      names(res) <- names(parameters)
    } else {
      rownames(mathessian) <- colnames(mathessian) <- names(x)
      rh <- SEfromHessian(mathessian, hessian=TRUE)
      res <- rh$SE
      result$hessian <- rh$hessian
    }
    
    if (any(is.na(res))) {
      warning("Problem in Hessian estimation. Confidence interval will not be available.")
      # res <- rep(NA, length(parameters))
    }
  }
  
  result$SE <- res 
  result$data <- temperatures
  
  # Avant *5. Correction du 17/7/2012
  k <- length(parameters)
  n <- unname(result$data$IndiceT["NbTS"])
  if ((all(!is.na(hatchling.metricuse[, "SD"]))) | (!is.na(result$par["SD"]))) {
    result$AIC <- 2*result$value+2*k
    result$AICc <- result$AIC + (2*k*(k+1))/(n-k-1)
    result$BIC <- 2*result$value + k * log(n)
  } else {
    result$AIC <- n*log(result$value/n)+2*k
    result$AICc <- result$AIC + (2*k*(k+1))/(n-k-1)
    result$BIC <- n*log(result$value/n) + k * log(n)
  }
  
  result$hatchling.metric <- hatchling.metricuse
  result$integral <- integral
  result$derivate <- derivate
  result$M0 <- M0
  # Je stocke aussi les paramtres fixe-16/7/2012
  result$fixed.parameters <- fixed.parameters
  # 29/1/2014
  result$weight <- weight
  
  
  result <- addS3Class(result, "NestsResult")
  
  return(result)
  
}
