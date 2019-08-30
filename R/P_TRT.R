#' P_TRT estimates the transitional range of temperatures based on a set of parameters
#' @title Estimate the transitional range of temperatures based on a set of parameters
#' @author Marc Girondot
#' @return A list with a P_TRT object containing a matrix with lower and higher bounds for TRT, TRT and P and a P_TRT_quantiles object with quantiles for each and a sexratio_quantiles object
#' @param x Set of parameters or a result of tsd()
#' @param resultmcmc A result of tsd_MHmcmc()
#' @param chain What chain to be used is resultmcmc is provided
#' @param temperatures If provided returns the sex ratio and its quantiles for each temperature
#' @param durations If provided returns the sex ratio and its quantiles for each duration
#' @param SD.temperatures SD of temperatures
#' @param SD.durations SD of durations
#' @param replicate.CI If a result of tsd() is provided, use replicate.CI replicates from the hessian matrix to estimate CI 
#' @param equation What equation should be used. Must be provided if x is not a result of tsd()
#' @param l Sex ratio limits to define TRT are l and 1-l. If NULL, TRT is not estimated.
#' @param probs Probabilities used to estimate quantiles
#' @param warn Do the warnings must be shown ? TRUE or FALSE
#' @description Estimate the parameters that best describe the thermal reaction norm for sex ratio when temperature-dependent sex determination occurs.\cr
#' It can be used also to evaluate the relationship between incubation duration and sex ratio.\cr
#' The parameter l was defined in Girondot (1999). The TRT is defined from the difference between the two boundary temperatures giving sex ratios of l and 1 - l.\cr
#' For logistic equation, exact value is used and precision iterations are used for other equations.\cr
#' In Girondot (1999), l was 0.05 and then the TRT was defined as being the range of temperatures producing from 5% to 95% of each sex.\cr
#' If l is null, TRT is not estimated and only sex ratio is estimated.\cr
#' @references Girondot, M. 1999. Statistical description of temperature-dependent sex determination using maximum likelihood. Evolutionary Ecology Research, 1, 479-486.
#' @references Godfrey, M.H., Delmas, V., Girondot, M., 2003. Assessment of patterns of temperature-dependent sex determination using maximum likelihood model selection. Ecoscience 10, 265-272.
#' @references Hulin, V., Delmas, V., Girondot, M., Godfrey, M.H., Guillon, J.-M., 2009. Temperature-dependent sex determination and global change: are some species at greater risk? Oecologia 160, 493-506.
#' @family Functions for temperature-dependent sex determination
#' @examples
#' \dontrun{
#' library("embryogrowth")
#' CC_AtlanticSW <- subset(DatabaseTSD, RMU=="Atlantic, SW" & 
#'                           Species=="Caretta caretta" & Sexed!=0)
#' tsdL <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="logistic"))
#' P_TRT(tsdL)
#' P_TRT(tsdL, replicate.CI=1000)
#' P_TRT(tsdL, replicate.CI=1000, temperatures=20:35)
#' P_TRT_out <- P_TRT(tsdL, replicate.CI=1000, temperatures=c(T1=20, T2=35))
#' attributes(P_TRT_out$sexratio_quantiles)$temperatures
#' P_TRT(tsdL$par, equation="logistic")
#' pMCMC <- tsd_MHmcmc_p(tsdL, accept=TRUE)
#' # Take care, it can be very long
#' result_mcmc_tsd <- tsd_MHmcmc(result=tsdL, 
#' 		parametersMCMC=pMCMC, n.iter=10000, n.chains = 1,  
#' 		n.adapt = 0, thin=1, trace=FALSE, adaptive=TRUE)
#' P_TRT(result_mcmc_tsd, equation="logistic")
#' }
#' @export

P_TRT <- function(x=NULL, resultmcmc=NULL, chain=1, equation=NULL, l=0.05, 
                  replicate.CI=NULL, temperatures=NULL, durations=NULL,
                  SD.temperatures= NULL, SD.durations=NULL,
                  probs=c(0.025, 0.5, 0.975), 
                  warn=TRUE) {
  
  modelTSD <- getFromNamespace(".modelTSD", ns="embryogrowth")
  
  # resultmcmc=NULL; chain=1; equation=NULL; l=0.05; replicate.CI=NULL; temperatures=NULL; durations=NULL; SD.temperatures= NULL; SD.durations=NULL; probs=c(0.025, 0.5, 0.975); warn=TRUE
  if (is.null(x) & is.null(resultmcmc)) stop("Or x or resultmcmc must be provided")
  if (!is.null(temperatures) & !is.null(durations)) stop("Temperatures and durations cannot be provided at the same time")
  
  if (class(x) == "tsd") {
    type <- x$type
  } else {
    if (!is.null(temperatures)) {
      type <- "temperature"
    } else {
      type <- "duration"
    }
  }
  
  if (!is.null(replicate.CI)) if (replicate.CI == 0) replicate.CI <- NULL
  
  temperatures <- c(temperatures, durations)
  SD <- c(SD.temperatures, SD.durations)
  
  if (is.null(SD) & !is.null(temperatures)) SD <- rep(0, length(temperatures))
  
  if (length(SD) != length(temperatures)) stop("Same number of SD and temperatures or durations must be provided")
  
  par <- NULL
  
  # TRT.limits <- c(1, 90)
  precision <- 15
  
  if (class(x)=="mcmcComposite" | class(resultmcmc)=="mcmcComposite") {
    if (class(x)=="mcmcComposite") {
      par <- x$resultMCMC[[chain]]
      if (is.null(equation)) equation <- x$equation
    }
    if (class(resultmcmc)=="mcmcComposite") {
      par <- resultmcmc$resultMCMC[[chain]]
      if (class(x)=="tsd") equation <- x$equation
      if (is.null(equation)) equation <- resultmcmc$equation
    }
  } else {
    if (class(x)=="tsd") {
      equation <- x$equation
      if (!is.null(replicate.CI) & (!is.null(x$hessian))) {
        if (requireNamespace("lmf")) {
          partot <- matrix(x$par, nrow=1)
          colnames(partot) <- names(x$par)
          partot <- partot[-1, ]
          
          errorsignS <- 0
          errorsignK1 <- 0
          errorsignK2 <- 0
          errorsignK <- 0
          errorsigncpt <- 0
          
          vcov <- solve(x$hessian)
          
          repeat {
            
            # 2019-05-31 : if replicate.CI == 1, renvoie quand même un nombre random
            par <- getFromNamespace("rmnorm", ns="lmf")(n = replicate.CI-nrow(partot), mean = x$par, vcov)
            # par <- getFromNamespace("rmnorm", ns="lmf")(n = replicate.CI-1, mean = x$par, vcov)
            # par <- rbind(x$par, par)
            if (!is.matrix(par)) {
              par <- matrix(par, nrow = 1)
              colnames(par) <- names(x$par)
            }
            
            errorsigncpt <- errorsigncpt + replicate.CI-nrow(partot)
            
            # Je retire quand S change de signe
            # Là je ne le retire pas mais je l'inverse !
            if (any(sign(par[, "S"]) != sign(x$par["S"]))) {
              errorsignS <- errorsignS + sum(sign(par[, "S"]) != sign(x$par["S"]))
              par <- par[sign(par[, "S"]) == sign(x$par["S"]), , drop = FALSE]
            }
            
            # Le K est forcément A-symmetrical, donc le pivot est 0
            if (!is.na(x$par["K"])) {
              if (nrow(par) != 0) {
                if (any(sign(par[, "K"]) != sign(x$par["K"]))) {
                  errorsignK <- errorsignK + sum(sign(par[, "K"]) != sign(x$par["K"]))
                  par <- par[sign(par[, "K"]) == sign(x$par["K"]), , drop = FALSE]
                }
              }
            }

            # K1 et K2 sont forcément de flexit, donc le pivot est 1
            if (!is.na(x$par["K1"])) {
              if (nrow(par) != 0) {
                if (any(sign(par[, "K1"] - 1) != sign(x$par["K1"] - 1))) {
                  errorsignK1 <- errorsignK1 + sum(sign(par[, "K1"] - 1) != sign(x$par["K1"] - 1))
                  par <- par[sign(par[, "K1"] - 1) == sign(x$par["K1"] - 1), , drop = FALSE]
                }
              }
            }
            
            if (!is.na(x$par["K2"])) {
              if (nrow(par) != 0) {
                if (any(sign(par[, "K2"] - 1) != sign(x$par["K2"] - 1))) {
                  errorsignK2 <- errorsignK2 + sum(sign(par[, "K2"] - 1) != sign(x$par["K2"] - 1))
                  par <- par[sign(par[, "K2"] - 1) == sign(x$par["K2"] - 1), , drop = FALSE]
                }
              }
            }
            
            
            if (nrow(par) != 0) {
              partot <- rbind(partot, par)
            }
            if (nrow(partot) == replicate.CI) break
          }
          
          if (warn) {
            if (errorsignK != 0) {
              message(paste("SE for K too high in", errorsignK, "cases out of", errorsigncpt))
            }
            if (errorsignS != 0) {
              message(paste("SE for S too high in", errorsignS, "cases out of", errorsigncpt))
            }
            if (errorsignK1 != 0) {
              message(paste("SE for K1 too high in", errorsignK1, "cases out of", errorsigncpt))
            }
            if (errorsignK2 != 0) {
              message(paste("SE for K2 too high in", errorsignK2, "cases out of", errorsigncpt))
            }
            if ((errorsignK + errorsignK2 + errorsignK1 + errorsignS) != 0) {
              warning("Use results with caution; it is probably better to use MCMC")
            }
          }
          par <- partot
          
          # colK <- which(colnames(par)=="K")
          # if (!identical(colK, integer(0))) {
          #   par <- par[sign(x$par[colK])==sign(par[, colK]), ]
          # }
        } else {
          warning("The package lmf is required to estimate confidence interval.")
          par <- t(as.matrix(x$par))
        }
      } else {
        par <- t(as.matrix(x$par))
      }
      
    } else {
      if (class(x)=="numeric") {
        par <- t(as.matrix(x))
      }
    }
  }
  
  if (is.null(par))  stop("I don't understand the format of x parameter")
  
  # print(equation)
  
  if (is.null(equation)) stop("equation parameter must be provided")
  equation <- tolower(equation)
  
  if (!is.null(temperatures)) {
    srT <- apply(X = par, MARGIN=1, function(xpar) {
      p <- modelTSD(xpar, rnorm(length(temperatures), temperatures, SD), equation)
      if (any(is.na(p))) {
        warning("Error for this set of parameters")
        print(d(xpar))
      }
      return(p)
    })
  } else {
    srT <- NULL
  }
  
  if (!is.null(l)) {
    
    # obj <- function(temperature, xpar, equation, objective) {
    #   p <- modelTSD(xpar, temperature, equation)
    #   # print(paste(temperature, unname(p), abs(p-objective)))
    #   return(abs(p-objective))
    # }
    
    ret <- apply(X = par, MARGIN=1, function(xpar) {
      # for (j in 1:nrow(par)) {
      #   print(j)
      #   xpar <- par[j, ]
      # print(xpar)
      outr <- NULL
      if (equation == "gsd") {
        limit.low.TRT <- -Inf
        limit.high.TRT <- +Inf
        outr <- c(lower.limit.TRT=unname(min(c(limit.low.TRT, limit.high.TRT))), 
                  higher.limit.TRT=unname(max(c(limit.low.TRT, limit.high.TRT))), 
                  TRT=unname(abs(limit.high.TRT-limit.low.TRT)), 
                  PT=unname(xpar["P"]))
      }
      if (equation=="logistic") {
        # sr <- 1/(1+exp(1/S)(P-T)) 
        # (1/l) - 1 <- exp(1/S)(P-T)
        # log((1/l)-1) <- (1/S)(P-T)
        # log((1/l)-1)*S <- P-T
        limit.low.TRT <- unname(xpar["P"]-log((1/l)-1)*xpar["S"])
        limit.high.TRT <- unname(xpar["P"]-log((1/(1-l))-1)*xpar["S"])
        outr <- c(lower.limit.TRT=unname(min(c(limit.low.TRT, limit.high.TRT))), 
                  higher.limit.TRT=unname(max(c(limit.low.TRT, limit.high.TRT))), 
                  TRT=unname(abs(limit.high.TRT-limit.low.TRT)), 
                  PT=unname(xpar["P"]))
      }
      if (equation=="flexit") {
        P <- xpar["P"]
        S <- xpar["S"]
        K1 <- xpar["K1"]
        K2 <- xpar["K2"]
        
        if (K1 == 0) {
          K1 <- 1e-9
        }
        if (K2 == 0) {
          K2 <- 1e-9
        }
        
        if (is.infinite(2^(K1))) {
          S1 <- K1*S
          K1 <- sign(K1)*500
        } else {
          S1 <- (2^(K1 - 1)*K1*S)/(2^(K1) - 1)
        }
        
        
        if (is.infinite(2^(K2))) {
          S2 <- K2*S
          K2 <- sign(K2)*500
        } else {
          S2 <- (2^(K2 - 1)*K2*S)/(2^(K2) - 1)
        }
        
        # (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)) = (1-l)
        # (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(1/K1)) = 1/(1-l)
        # (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x))) = (1/(1-l)) ^ K1
        # (2^K1 - 1) *  exp(4 * S1 * (P - x)) = (1/(1-l)) ^ K1 - 1
        # exp(4 * S1 * (P - x)) = ((1/(1-l)) ^ K1 - 1)/(2^K1 - 1)
        # 4 * S1 * (P - x) = log(((1/(1-l)) ^ K1 - 1)/(2^K1 - 1))
        # P - x = log(((1/(1-l)) ^ K1 - 1)/(2^K1 - 1))/(4 * S1)
        # x= P -  log(((1/(1-l)) ^ K1 - 1)/(2^K1 - 1))/(4 * S1)
        
        limit.low.TRT <- P -  log(((1/(1-l)) ^ K1 - 1)/(2^K1 - 1))/(4 * S1)
        
        # 1-(1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2) = l
        # (1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2) = (1-l)
        # (1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(1/K2) = 1/(1-l)
        # 1 + (2^K2 - 1) * exp(4 * S2 * (x - P)) = (1/(1-l))^K2
        # (2^K2 - 1) * exp(4 * S2 * (x - P)) = (1/(1-l))^K2 - 1
        # exp(4 * S2 * (x - P)) = ((1/(1-l))^K2 - 1)/(2^K2 - 1)
        # 4 * S2 * (x - P)) = log(((1/(1-l))^K2 - 1)/(2^K2 - 1))
        # x - P = 1/(4 * S2) * log(((1/(1-l))^K2 - 1)/(2^K2 - 1))
        # x = 1/(4 * S2) * log(((1/(1-l))^K2 - 1)/(2^K2 - 1)) + P
        
        limit.high.TRT <- 1/(4 * S2) * log(((1/(1-l))^K2 - 1)/(2^K2 - 1)) + P
        
        outr <- c(lower.limit.TRT=unname(limit.low.TRT), 
                  higher.limit.TRT=unname(limit.high.TRT), 
                  TRT=unname(limit.high.TRT-limit.low.TRT), 
                  PT=unname(P))
      }
      
      if (equation=="a-logistic") {
        P <- xpar["P"]
        S <- xpar["S"]
        K <- xpar["K"]

        # (1 + (2^exp(K) - 1) * exp((1/S) * (P - x)))^(-1/exp(K)) =  1 - l
        # (1 + (2^exp(K) - 1) * exp((1/S) * (P - x)))^(1/exp(K)) =  1/(1 - l)
        # 1 + (2^exp(K) - 1) * exp((1/S) * (P - x)) =  (1/(1 - l))^exp(K)
        # (2^exp(K) - 1) * exp((1/S) * (P - x)) =  ((1/(1 - l))^exp(K) - 1)
        # exp((1/S) * (P - x)) =  ((1/(1 - l))^exp(K) - 1) / (2^exp(K) - 1)
        # (1/S) * (P - x) = log( ((1/(1 - l))^exp(K) - 1) / (2^exp(K) - 1) )
        # P - x = S * log( ((1/(1 - l))^exp(K) - 1) / (2^exp(K) - 1) )
        # x = P - S * log( ((1/(1 - l))^exp(K) - 1) / (2^exp(K) - 1) )
        
        
        limit.low.TRT <- P - S * log( ((1/(1 - l))^ exp(K) - 1) / (2^exp(K) - 1) )
        
        if (is.infinite(limit.low.TRT) & (K > 0)) {
          # exp((-1/exp(K))*(exp(K)*log(2)+(1/S)*(P-x))) =  1 - l
          # (-1/exp(K))*(exp(K)*log(2)+(1/S)*(P-x)) =  log(1 - l)
          # exp(K)*log(2)+(1/S)*(P-x) =  log(1 - l)/(-1/exp(K))
          # log(2)+(1/S)*(P-x) =  (log(1 - l)/(-1/exp(K)))/exp(K)
          # (1/S)*(P-x) =  ((log(1 - l)/(-1/exp(K)))/exp(K))-log(2)
          # P-x =  (((log(1 - l)/(-1/exp(K)))/exp(K))-log(2))*S
          # x =  P-(((log(1 - l)/(-1/exp(K)))/exp(K))-log(2))*S
          
          limit.low.TRT <- P-(((log(1 - l)/(-1/exp(K)))/exp(K))-log(2))*S
        }
        
        if ((!is.finite(limit.low.TRT)) & (K < 0)) {
          # x = P - S * log( ((1/(1 - l))^exp(K) - 1) / (2^exp(K) - 1) )
          # Kl = (1/(1 - l))
          # eK = exp(K)
          # x = P - S * log( (Kl^eK - 1) / (2^eK - 1) )
          # (Kl^eK - 1) / (2^eK - 1)
          # Approximation linéaire
          # Dérivée de Kl^eK - 1 = eK * Kl^(eK -1) * eK = eK^2 Kl^(eK -1)
          # Dérivée de 2^eK - 1 = eK * 2^(eK -1) * eK = eK^2 2^(eK -1)
          # (Kl^eK - 1) / (2^eK - 1) = Kl^(eK -1) / 2^(eK -1)
          # (Kl / 2)^(eK - 1)
          # x = P - S * log( (Kl / 2)^(eK - 1) )
          # x = P - S * (eK - 1)* log( (Kl / 2))
          # x = P - S * (exp(K) - 1)* log( ((1/(1 - l))/ 2))
          
          limit.low.TRT <- P - S * (exp(K) - 1)* log( ((1/(1 - l))/ 2))
        }
        
        # (1 + (2^exp(K) - 1) * exp((1/S) * (P - x)))^(-1/exp(K)) =  l
        # (1 + (2^exp(K) - 1) * exp((1/S) * (P - x)))^(1/exp(K)) =  1/l
        # 1 + (2^exp(K) - 1) * exp((1/S) * (P - x)) =  (1/l)^exp(K)
        # (2^exp(K) - 1) * exp((1/S) * (P - x)) =  ((1/l)^exp(K) - 1)
        # exp((1/S) * (P - x)) =  ((1/l)^exp(K) - 1) / (2^exp(K) - 1)
        # (1/S) * (P - x) = log( ((1/l)^exp(K) - 1) / (2^exp(K) - 1) )
        # P - x = S * log( ((1/l)^exp(K) - 1) / (2^exp(K) - 1) )
        # x = P - S * log( ((1/l)^exp(K) - 1) / (2^exp(K) - 1) )
        
        
        limit.high.TRT <- P - S * log( ((1 / l)^ exp(K) - 1) / (2^exp(K) - 1) )
        
        if (is.infinite(limit.high.TRT) & (K > 0)) {
          # exp((-1/exp(K))*(exp(K)*log(2)+(1/S)*(P-x))) =  l
          # (-1/exp(K))*(exp(K)*log(2)+(1/S)*(P-x)) =  log(l)
          # exp(K)*log(2)+(1/S)*(P-x) =  log(l)/(-1/exp(K))
          # log(2)+(1/S)*(P-x) =  (log(l)/(-1/exp(K)))/exp(K)
          # (1/S)*(P-x) =  ((log(l)/(-1/exp(K)))/exp(K))-log(2)
          # P-x =  (((log(l)/(-1/exp(K)))/exp(K))-log(2))*S
          # x =  P-(((log(l)/(-1/exp(K)))/exp(K))-log(2))*S
          
          limit.high.TRT <- P-(((log(l)/(-1/exp(K)))/exp(K))-log(2))*S
        }
        
        if ((!is.finite(limit.high.TRT)) & (K < 0)) {
          # x = P - S * log( ((1/(1 - l))^exp(K) - 1) / (2^exp(K) - 1) )
          # Kl = (1/(1 - l))
          # eK = exp(K)
          # x = P - S * log( (Kl^eK - 1) / (2^eK - 1) )
          # (Kl^eK - 1) / (2^eK - 1)
          # Approximation linéaire
          # Dérivée de Kl^eK - 1 = eK * Kl^(eK -1) * eK = eK^2 Kl^(eK -1)
          # Dérivée de 2^eK - 1 = eK * 2^(eK -1) * eK = eK^2 2^(eK -1)
          # (Kl^eK - 1) / (2^eK - 1) = Kl^(eK -1) / 2^(eK -1)
          # (Kl / 2)^(eK - 1)
          # x = P - S * log( (Kl / 2)^(eK - 1) )
          # x = P - S * (eK - 1)* log( (Kl / 2))
          # x = P - S * (exp(K) - 1)* log( ((1/(1 - l))/ 2))
          
          limit.high.TRT <- P - S * (exp(K) - 1)* log( ((1/ l)/ 2))
        }
        
        if (is.infinite(limit.high.TRT) | is.infinite(limit.low.TRT)) {
          print(d(c(P=P, S=S, K=K)))
        }
        
        outr <- c(lower.limit.TRT=unname(limit.low.TRT), 
                  higher.limit.TRT=unname(limit.high.TRT), 
                  TRT=unname(limit.high.TRT-limit.low.TRT), 
                  PT=unname(P))
      }
      
      if (equation=="hill") {
        P <- xpar["P"]
        S <- xpar["S"]
       
        # 1/(1+exp((1/par["S"])*(log(par["P"])-log(temperatures)))) = 1 - l
        # 1+exp((1/par["S"])*(log(par["P"])-log(temperatures))) = 1 / (1 - l)
        # exp((1/par["S"])*(log(par["P"])-log(temperatures))) = (1 / (1 - l))-1
        # (1/par["S"])*(log(par["P"])-log(temperatures)) = log((1 / (1 - l))-1)
        # log(par["P"])-log(temperatures) = log((1 / (1 - l))-1) * par["S"]
        # log(par["P"])-log(temperatures) = log((1 / (1 - l))-1) * par["S"]
        # log(temperatures) = log(par["P"]) - (log((1 / (1 - l))-1) * par["S"])
        # temperatures = exp(log(par["P"]) - (log((1 / (1 - l))-1) * par["S"]))
        
        limit.low.TRT <- exp(log(P) - (log((1 / (1 - l))-1) * S))
        limit.high.TRT <- exp(log(P) - (log((1 / l)-1) * S))
        
        outr <- c(lower.limit.TRT=unname(limit.low.TRT), 
                  higher.limit.TRT=unname(limit.high.TRT), 
                  TRT=unname(limit.high.TRT-limit.low.TRT), 
                  PT=unname(P))
      }
      
      # Les autres modèles
      if (is.null(outr)) {
        outr <- c(lower.limit.TRT=NA, 
                  higher.limit.TRT=NA, 
                  TRT=NA, 
                  PT=NA)
      }
      return(outr)
    }
    )
    
    if (any(is.na(ret)) & warn) warning("Something strange occurs; I cannot estimate TRT limits for some replicates")
    
    ret <- t(ret)
    pr <- apply(ret, MARGIN=2, function(xxx) quantile(xxx, probs = probs, na.rm=TRUE))
  } else {
    ret <- NULL
    pr <- NULL
  }
  
  if (!is.null(srT)) {
    # Là c'est une température
    if (is.null(dim(srT))) {
      prsrT <- quantile(srT, probs = probs, na.rm = TRUE)
      prsrT <- as.matrix(prsrT)
    } else {
      prsrT <- apply(t(srT), MARGIN = 2, function(xxx) quantile(xxx, probs = probs, na.rm = TRUE))
      if (class(prsrT) == "numeric") {
        prsrT <- matrix(prsrT, nrow = 1)
        rownames(prsrT) <- paste0(as.character(probs*100), "%")
      }
    }
    if (!is.null(names(temperatures))) {
      colnames(prsrT) <- names(temperatures)
    } else {
      colnames(prsrT) <- as.character(temperatures)
    }
    if (type =="duration") {
      attributes(prsrT) <- c(attributes(prsrT), durations=list(unname(temperatures)))
    } else {
      attributes(prsrT) <- c(attributes(prsrT), temperatures=list(unname(temperatures)))
    }
  } else {
    prsrT <- NULL
  }
  return(list(P_TRT=ret, P_TRT_quantiles=pr, sexratio_quantiles=prsrT))
}
