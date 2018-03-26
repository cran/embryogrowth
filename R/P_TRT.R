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
                warn=FALSE) {
  
  modelTSD <- getFromNamespace(".modelTSD", ns="embryogrowth")
  
  # resultmcmc=NULL; chain=1; equation=NULL; l=0.05; replicate.CI=NULL; temperatures=NULL; durations=NULL; SD.temperatures= NULL; SD.durations=NULL; probs=c(0.025, 0.5, 0.975); warn=FALSE
  if (is.null(x) & is.null(resultmcmc)) stop("Or x or resultmcmc must be provided")
  if (!is.null(temperatures) & !is.null(durations)) stop("Temperatures and durations cannot be provided at the same time")
  
  if (class(x)=="tsd") {
    type <- x$type
    # if (is.null(temperatures) & is.null(durations)) {
    #   temperatures <- x$temperatures
    # }
  } else {
    if (!is.null(temperatures)) {
      type <- "temperature"
    } else {
      type <- "duration"
    }
  }
  
  if (!is.null(replicate.CI)) if (replicate.CI <= 1) replicate.CI <- NULL
  
  temperatures <- c(temperatures, durations)
  SD <- c(SD.temperatures, SD.durations)
  
  if (is.null(SD) & !is.null(temperatures)) SD <- rep(0, length(temperatures))
  
  if (length(SD) != length(temperatures)) stop("Same number of SD and temperatures or durations must be provided")
  
  par <- NULL
  
  # TRT.limits <- c(1, 90)
  precision <- 15
  
  if (class(x)=="mcmcComposite" | class(resultmcmc)=="mcmcComposite") {
    if (class(x)=="mcmcComposite") par <- x$resultMCMC[[chain]]
    if (class(resultmcmc)=="mcmcComposite") {
      par <- resultmcmc$resultMCMC[[chain]]
    if (class(x)=="tsd") equation <- x$equation
    }
  } else {
  if (class(x)=="tsd") {
    equation <- x$equation
    if (!is.null(replicate.CI) & (!is.null(x$hessian))) {
      if (requireNamespace("lmf")) {
      vcov <- solve(x$hessian)
      par <- getFromNamespace("rmnorm", ns="lmf")(n = replicate.CI-1, mean = x$par, vcov)
      par <- rbind(x$par, par)
      colnames(par) <- names(x$par)
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
      p <- modelTSD(xpar, rnorm(length(temperatures)+3, c(xpar["P"]-20, xpar["P"], xpar["P"]+20, temperatures), c(0, 0, 0, SD)), equation)

      if (all(is.finite(p))) {
        # C'est quoi ça ?
      if ((abs(p[1]-p[2])<0.1) | (abs(p[2] - p[3])<0.1) | (abs(p[2]-0.5) > 0.1)) { #| (p[1] < 0.95) | (p[3] > 0.05)
        return(rep(NA, length(temperatures)))
      } else {
        return(p[-(1:3)])
      }
      } else {
        return(rep(NA, length(temperatures)))
      }
    })
  } else {
    srT <- NULL
  }
  
  if (!is.null(l)) {
  
 ret <- apply(X = par, MARGIN=1, function(xpar) {
     # for (j in 1:nrow(par)) {
     #   print(j)
     #   xpar <- par[j, ]
  if (equation=="gsd") {
    limit.low.TRT <- -Inf
    limit.high.TRT <- +Inf
    outr <- c(lower.limit.TRT=unname(min(c(limit.low.TRT, limit.high.TRT))), 
              higher.limit.TRT=unname(max(c(limit.low.TRT, limit.high.TRT))), 
              TRT=unname(abs(limit.high.TRT-limit.low.TRT)), 
              PT=unname(xpar["P"]))
    
  } else {
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
    
  
    } else {

    limit.low.TRT <- c(xpar["P"]-20, xpar["P"], xpar["P"]+20)
    limit.high.TRT <- c(xpar["P"]-20, xpar["P"], xpar["P"]+20)
    p <- modelTSD(xpar, limit.low.TRT, equation)
    if (all(is.finite(p))) {
      if ((abs(p[1]-p[2])<0.1) | (abs(p[2] - p[3])<0.1) | (abs(p[2]-0.5) > 0.1)) {
      outr <- c(lower.limit.TRT=NA, 
                higher.limit.TRT=NA, 
                TRT=NA, 
                PT=NA)
    } else {
      cpt <- 1
      
      repeat {
        p <- modelTSD(xpar, limit.low.TRT, equation)
        if (p[2]<(1-l)) {
          limit.low.TRT <- c(limit.low.TRT[1], limit.low.TRT[1]+(limit.low.TRT[2]-limit.low.TRT[1])/2, limit.low.TRT[2])
        } else {
          limit.low.TRT<- c(limit.low.TRT[2], limit.low.TRT[2]+(limit.low.TRT[3]-limit.low.TRT[2])/2, limit.low.TRT[3])
        }
        p <- modelTSD(xpar, limit.high.TRT, equation)
        if (p[2]<l) {
          limit.high.TRT <- c(limit.high.TRT[1], limit.high.TRT[1]+(limit.high.TRT[2]-limit.high.TRT[1])/2, limit.high.TRT[2])
        } else {
          limit.high.TRT <- c(limit.high.TRT[2], limit.high.TRT[2]+(limit.high.TRT[3]-limit.high.TRT[2])/2, limit.high.TRT[3])
        }
        
        cpt <- cpt + 1
        if ((cpt>=15) | ((abs(limit.high.TRT[2]-limit.high.TRT[1])<1E-3) & (abs(limit.low.TRT[2]-limit.low.TRT[1])<1E-3))) break
      }
    
      limit.low.TRT <- limit.low.TRT[2]
      limit.high.TRT <- limit.high.TRT[2]
    
  
    outr <- c(lower.limit.TRT=unname(min(c(limit.low.TRT, limit.high.TRT))), 
              higher.limit.TRT=unname(max(c(limit.low.TRT, limit.high.TRT))), 
              TRT=unname(abs(limit.high.TRT-limit.low.TRT)), 
              PT=unname(xpar["P"]))
    }
    } else {
      outr <- c(lower.limit.TRT=NA, 
                higher.limit.TRT=NA, 
                TRT=NA, 
                PT=NA)
    }
  }
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
