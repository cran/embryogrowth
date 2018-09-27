#' uncertainty.datalogger Calculate the uncertainty of the average temperature calculated using data gathered by a data logger.
#' @title Uncertainty of average temperatures obtained using temperature data logger
#' @author Marc Girondot
#' @return The function will return the uncertainty of the average temperature for the considered period as being the 95% range where the true average temperature should be.
#' @param max.time being the maximum time to record in minutes
#' @param sample.rate The sample rates in minutes
#' @param accuracy The accuracy of the data logger in °C
#' @param resolution The resolution of the data logger in °C
#' @param replicates The number of replicates to estimate uncertainty.
#' @param method The fonction that will be used to return the uncertainty.
#' @description Calculate the uncertainty of average temperature dependent on the 
#' characteristics of a data logger and sampling rate.\cr
#' The temperature is supposed to be uniformaly distributed with min and max 
#' being -accuracy and +accuracy.
#' @family Data loggers utilities
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' # Exemple using the hypothesis of Gaussian distribution
#' uncertainty.datalogger(sample.rate=30, accuracy=1, resolution=0.5, 
#'                 method=function(x) {2*qnorm(0.975)*sd(x)})
#' # Example without hypothesis about distribution, using quantiles
#' uncertainty.datalogger(sample.rate=30, accuracy=1, resolution=0.5, 
#'                 method=function(x) {quantile(x, probs=c(0.975))-
#'                                               quantile(x, probs=c(0.025))})
#' par(mar=c(4, 4, 1, 1))
#' plot(x=10:120, uncertainty.datalogger(sample.rate=10:120, 
#'                                       accuracy=0.5, 
#'                                       resolution=1), 
#'      las=1, bty="n", type="l", 
#'      xlab="Sample rate in minutes", 
#'      ylab=expression("Uncertainty in "*degree*"C"), 
#'      ylim=c(0, 0.15), xlim=c(0, 120))  
#' lines(x=10:120, uncertainty.datalogger(sample.rate=10:120, 
#'                                             accuracy=1, 
#'                                             resolution=0.5), col="red")
#' lines(x=10:120, uncertainty.datalogger(sample.rate=10:120, 
#'                                        accuracy=1, 
#'                                        resolution=1), col="blue")
#' lines(x=10:120, uncertainty.datalogger(sample.rate=10:120, 
#'                                        accuracy=0.5, 
#'                                        resolution=0.5), col="yellow")
#' legend("topleft", legend=c("Accuracy=0.5, resolution=0.5", 
#'                            "Accuracy=0.5, resolution=1", 
#'                            "Accuracy=1, resolution=0.5", 
#'                            "Accuracy=1, resolution=1"), lty=1, 
#'        col=c("yellow", "black", "red", "blue"), 
#'        cex=0.6)
#' }
#' @export


uncertainty.datalogger <- function(max.time=0, 
                sample.rate=0, 
                accuracy=0.5, 
                resolution=1, 
                replicates=10000, 
                method=function(x) {2*qnorm(0.975)*sd(x)}) {
  
  # max.time=10*24*60 
  # sample.rate=15
  # accuracy=0.5 
  # resolution=1
  # replicates=1000
  # method=function(x) {2*qnorm(0.975)*sd(x)}
  
  if (replicates < 2) {
    stop("replicates must be at least 2.")
  }
  
  GT <- NULL
  
  for (sr in sample.rate) {
    
  time <- seq(from=0, to=max.time, by=sr)
  sampledT <- NULL
  re <- seq(from=-1, to=1, length.out = replicates)
  for (i in 1:replicates) {
    temp <- runif(1, min=25, max=30)+time*re[i]*0.0002
    temp1 <- temp+runif(length(temp), min=-accuracy, max=accuracy)
    temp2 <- floor((temp1+resolution/2)*(1/resolution))*resolution
    sampledT <- c(sampledT, mean(temp2)-mean(temp))
  }
  GT <- c(GT, method(sampledT))
  }
  names(GT) <- as.character(sample.rate)
  
  return(GT)
}

