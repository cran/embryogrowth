#' ChangeSSM convert different forms of thermal norm of reaction
#' @title Generate set of parameters for different forms of thermal norm of reaction
#' @author Marc Girondot
#' @return A vector with parameters or a result object formatted with new parameters is result is non null
#' @param result A result obtained by searchR()
#' @param resultmcmc A result obtained by GRTRN_MHmcmc()
#' @param temperatures A vector with incubation temperatures in degrees Celsius
#' @param parameters A vector of parameters for model to be converted. Not necessary if result is provided.
#' @param initial.parameters NULL or a vector of parameters for initial model model to be fited
#' @param fixed.parameters NULL of a vector of parameters to be used but fixed
#' @param outmcmc What statistic will be estimated if a mcmc is provided. Can be "mean-sd" or "quantiles".
#' @param progressbar If TRUE, a progressbar is shown
#' @param ... A control list to be used with optim, see ?optim
#' @description Generate a set of parameters for thermal reaction norm model.\cr
#' If initial.parameters is NULL and resultmcmc is not NULL, it will generate parameters and SE based on the average of the curves.
#' @examples
#' \dontrun{
#' data(resultNest_6p_SSM)
#' x1 <- resultNest_6p_SSM$par
#' data(resultNest_4p_SSM)
#' x2 <- resultNest_4p_SSM$par
#' temperaturesC <- (200:350)/10
#' s <- ChangeSSM(temperatures=temperaturesC, parameters=x1, initial.parameters=x2)
#' sY <- plotR(resultNest_6p_SSM, ylim=c(0,3), col="black", curve = "ML")
#' plotR(resultNest_4p_SSM, col="red", scaleY=sY, new=FALSE)
#' plotR(s$par, col="green", scaleY=sY, new=FALSE, curve = "ML")
#' legend("topleft", legend=c("r function to mimic", "Initial new r function", 
#' "Fitted new r function"), lty=c(1, 1, 1), col=c("black", "red", "green"))
#' # Other example to fit anchored parameters
#' data(resultNest_4p_SSM)
#' x0 <- resultNest_4p_SSM$par
#' t <- hist(resultNest_4p_SSM, plot=FALSE)
#' x <- c(3.4, 3.6, 5.4, 5.6, 7.6, 7.5, 3.2)
#' names(x) <- seq(from=range(t$temperatures)[1], to=range(t$temperatures)[2], 
#'      length.out=7)
#' newx <- ChangeSSM(temperatures = (200:350)/10, parameters = x0, 
#'        initial.parameters = x, 
#'        control=list(maxit=5000))
#'  # Example on how to generate a set of SSM parameters from anchored parameters
#'  xanchor <- GenerateAnchor(nests=resultNest_4p_SSM)
#'  x <- resultNest_4p_SSM$par
#'  xanchor["294"] <- 0
#'  xanchor["308"] <- 2.3291035
#'  x <- ChangeSSM(parameters = xanchor,
#'                      initial.parameters = x, control=list(maxit=5000))
#'  sY <- plotR(resultNest_4p_SSM$par, ylim = c(0,3), curve="ML")
#'  plotR(xprime$par, col="red", scaleY=sY, new=FALSE, curve="ML") 
#'  legend("topleft", legend=c("Fitted parameters", "Constrainted parameters"), lty=1, 
#'         col=c("black", "red"))
#'  # Weibull model
#'  x <- ChangeSSM(temperatures = (200:350)/10,
#'                 parameters = resultNest_4p_SSM$par,
#'                 initial.parameters = structure(c(73, 300, 26), 
#'                                                .Names = c("k", "lambda", "scale")), 
#'                 control=list(maxit=1000))
#'  # normal asymmetric model
#'  x <- ChangeSSM(temperatures = (200:350)/10,
#'                parameters = resultNest_4p_SSM$par,
#'                initial.parameters = structure(c(3, 10, 8, 32), 
#'                .Names = c("Scale", "sdL", "sdH", "Peak")), 
#'                control=list(maxit=1000))
#'  # trigonometric model
#'  x <- ChangeSSM(temperatures = (200:350)/10,
#'                parameters = resultNest_4p_SSM$par,
#'                initial.parameters = structure(c(3, 20, 40, 32), 
#'                .Names = c("Max", "LengthB", "LengthE", "Peak")), 
#'                control=list(maxit=1000))
#'
#'  # example with a mcmc object, CI being 2.SD
#'  # Note the symmetric CI
#' data(resultNest_mcmc_4p_SSM)
#' new_result <- ChangeSSM(resultmcmc = resultNest_mcmc_4p_SSM, result = resultNest_4p_SSM,
#'                         temperatures = seq(from = 20, to = 35, by = 0.1), 
#'                         outmcmc = "mean-sd", 
#'                         initial.parameters = NULL)
#' 
#' plotR(new_result, ylim=c(0, 3), curve="ML")
#'  # example with a mcmc object, CI being defined by 2.5%-97.5% quantiles
#'  # Note the asymmetric CI
#' data(resultNest_mcmc_4p_SSM)
#' new_result <- ChangeSSM(resultmcmc = resultNest_mcmc_4p_SSM, result = resultNest_4p_SSM,
#'                         temperatures = seq(from = 20, to = 35, by = 0.1), 
#'                         outmcmc = "quantiles", 
#'                         initial.parameters = NULL)
#'  
#' plotR(new_result, ylim=c(0, 3), curve="ML")
#' plotR(new_result, ylim=c(0, 3), curve="ML quantiles")
#' 
#' # A little trick
#' # to convert SSM4 to SSM6, you can use:
#' 
#' x4 <- c('DHA' = 69.718935117894063, 
#'          'DHH' = 497.81709040501079, 
#'          'T12H' = 308.95543713889509, 
#'          'Rho25' = 255.24186073771696)
#' 
#' x6 <- c(x4["DHA"], 
#'.        x4["DHH"], 
#'.        'DHL' = x4[["DHH"]], 
#'.        'DT' = 0.5, 
#'.        'T12L' = x4[["T12H"]], 
#'.        x4['Rho25'])
#' }
#' @export

ChangeSSM <- function(result = NULL, 
                      resultmcmc = NULL,
                      temperatures = seq(from = 20, to = 35, by = 0.1), 
                      parameters = NULL, 
                      initial.parameters = NULL, 
                      fixed.parameters = NULL,
                      outmcmc = "quantiles", 
                      progressbar = TRUE, 
                      ...) {
  
  # result = NULL; resultmcmc = NULL; temperatures = seq(from = 20, to = 35, by = 0.1); parameters = NULL; fixed.parameters = NULL; initial.parameters = NULL; outmcmc = "quantiles"; progressbar = TRUE

    
  if (is.null(parameters) & is.null(result) & is.null(resultmcmc)) stop("Or parameters or result or resultmcmc must be supplied.")

  if (is.null(parameters)) {
    if (!is.null(result)) {
      parameters <- c(result$par, result$fixed.parameters)
    } else {
      parameters <- structure(as.numeric(resultmcmc$parametersMCMC$parameters[, "Init"]), .Names=rownames(resultmcmc$parametersMCMC$parameters))
    }
  }
  
  if (is.null(initial.parameters)) {
    
    xxw <- matrix(data = as.numeric(), ncol = length(temperatures), 
                  nrow = nrow(resultmcmc$resultMCMC[[1]]))
    if (progressbar) pb<-txtProgressBar(min=1, max=nrow(resultmcmc$resultMCMC[[1]]), style=3)
    
    for (i in 1:nrow(resultmcmc$resultMCMC[[1]])) {
      if (progressbar) 	setTxtProgressBar(pb, i)
      xxw[i, ] <- getFromNamespace(".SSM", ns="embryogrowth")(T=temperatures, 
                                                       parms=resultmcmc$resultMCMC[[1]][i,])[[1]]*1E5
    }
    
    rxxw <- apply(xxw, MARGIN = 2, 
                  FUN = function(x) c(mean(x), sd(x), quantile(x = x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)))
    
    if (outmcmc == "quantiles") {
      newp <- structure(rxxw[4, ], .Names = as.character(x = temperatures+273.15))
      newse <- matrix(data = c(rxxw[3, ], rxxw[5, ]), nrow=2, byrow = TRUE, 
                      dimnames = list(c("2.5%", "97.5%"),  as.character(x = temperatures+273.15)))
    } else {
      newp <- structure(rxxw[1, ], .Names = as.character(x = temperatures+273.15))
      newse <- structure(rxxw[2, ], .Names = as.character(x = temperatures+273.15))
    }
    
    if (!is.null(result)) {
      result$par <- newp
      result$SE <- newse
      return(result)
    } else {
      return(list(par=newp, SE=newse))
    }
    
  } else {
  
growth.rate <- getFromNamespace(".SSM", ns="embryogrowth")(273.15+temperatures, parameters)[[1]]*1E5

c <- list(...)

if (length(c)==0) {
  s <- optim(par=initial.parameters, fn=getFromNamespace(".fitSSM", ns="embryogrowth"), temperatures=temperatures, growth.rate=growth.rate, fixed.parameters=fixed.parameters)
} else {
  s <- optim(par=initial.parameters, fn=getFromNamespace(".fitSSM", ns="embryogrowth"), temperatures=temperatures, growth.rate=growth.rate, fixed.parameters=fixed.parameters, control=c[[1]])
}
  if (!is.null(result)) {
    result$par <- s$par
    return(result)
  } else {
    return(s)
  }
}

}

