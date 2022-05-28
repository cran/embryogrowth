#' predict.HatchingSuccess returns prediction based on a model fitted with HatchingSuccessfit()
#' @title Return prediction based on a model fitted with HatchingSuccess.fit()
#' @author Marc Girondot
#' @return Return a matrix with prediction based on a model fitted with HatchingSuccess.fit()
#' @param temperature A vector of temperatures.
#' @param probs Quantiles.
#' @param replicates Number of replicates to estimate the confidence interval.
#' @param object The return of a fit done with HatchingSuccess.fit().
#' @param resultmcmc Results obtained using HatchingSuccess.MHmcmc()
#' @param chain Chain to use in resultmcmc
#' @param ... Not used
#' @description Set of functions to study the hatching success.\cr
#' If replicates is 0, it returns only the fitted model.\cr
#' If replicates is null and resultmcmc is not null, it will use all the mcmc data.\cr
#' if replicates is lower than the number of iterations in resultmcmc, it will use sequence of data regularly thined.
#' @family Hatching success
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' totalIncubation_Cc <- subset(DatabaseTSD, 
#'                              Species=="Caretta caretta" & 
#'                                Note != "Sinusoidal pattern" & 
#'                                !is.na(Total) & Total != 0)
#' 
#' par <- c(S.low=0.5, S.high=0.3, 
#'          P.low=25, deltaP=10, MaxHS=0.8)
#'          
#' HatchingSuccess.lnL(x=par, data=totalIncubation_Cc)
#' 
#' g <- HatchingSuccess.fit(par=par, data=totalIncubation_Cc)
#' 
#' HatchingSuccess.lnL(par=g$par, data=totalIncubation_Cc)
#' 
#' plot(g)
#' }
#' @method predict HatchingSuccess
#' @export




predict.HatchingSuccess <- function(object, ..., 
                                    temperature=NULL, 
                                    probs=c(0.025, 0.5, 0.975), 
                                    replicates=NULL, resultmcmc=NULL, chain=1) {
  
  par <- c(object$par, object$fixed.parameters)
  
  # mathessian <- object$hessian
  if (is.null(replicates) & is.null(resultmcmc)) {
    replicates <- 0
  }
  
  if (is.null(temperature)) temperature <- object$data[, object$column.Incubation.temperature]
  
  if (!is.null(resultmcmc)) {
    
    par <- resultmcmc$resultMCMC[[chain]]
    if (!is.null(replicates)) {
      if (replicates<nrow(par)) {
        par <- par[round(seq(from=1, to=nrow(par), length.out =  replicates)), ]
        
      }
    }
    replicates <- nrow(par)
    CI <- matrix(data = NA , ncol=replicates, nrow=length(temperature))
    
    for (c in 1:replicates) {
      CI[, c] <- HatchingSuccess.model(par=par[c, ], temperature)
    }
    
  } else {
    
    if (replicates >0) {
      
      CI <- matrix(data = NA , ncol=replicates, nrow=length(temperature))
      
      if (requireNamespace("lmf")) {
        vcov <- solve(object$hessian)
        # 2019-05-31 : if replicate.CI == 1, renvoie quand même un nombre random
        par <- getFromNamespace("rmnorm", ns="lmf")(n = replicates, mean = object$par, vcov)
        # par <- getFromNamespace("rmnorm", ns="lmf")(n = replicate.CI-1, mean = x$par, vcov)
        # par <- rbind(x$par, par)
        if (!is.matrix(par)) {
          par <- matrix(par, nrow = 1)
          colnames(par) <- names(object$par)
        }
        
        for (c in 1:replicates) {
          parec <- par[c, ]
          names(parec) <- colnames(object$hessian)
          CI[, c] <- HatchingSuccess.model(par=parec, temperature)
        }
      } else {
        warning("The package lmf should be present to better estimate confidence interval taking into account covariances.")
        for (c in 1:replicates) {
          SE <- object$SE
          x <- structure(rnorm(n = 5, mean=par, sd=SE), .Names=names(par))
          CI[, c] <- HatchingSuccess.model(par=x, temperature)
        }
        
      }
    } else {
      CI <- matrix(data = NA , ncol=1, nrow=length(temperature))
      CI[, 1] <- HatchingSuccess.model(par, temperature)
    }
  }
  
  CIq <- apply(CI, MARGIN = 1, FUN = function(x) {quantile(x, probs = probs)})
  if (inherits(CIq, "numeric")) CIq <- matrix(data=CIq, ncol=length(temperature))
  colnames(CIq) <- as.character(temperature)
  
  return(CIq)
}

