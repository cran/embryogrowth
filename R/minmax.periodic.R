#' minmax.periodic search for minimum and maximum temperatures in periodic timeseries
#' @title Search for minimum and maximum temperatures in periodic timeseries
#' @author Marc Girondot
#' @return A data.frame with a column time, a column temperature and a column sd
#' @param time.minmax.daily A named vector with Min and Max being the time in the day with minimum and maximum temperatures
#' @param time.obs	A vector with the time at which temperatures are recorded
#' @param temp.obs	A vector with the recorded temperatures recorded
#' @param period The unit of day period (24 for hours, 24*60 for minutes)
#' @description Search for minimum and maximum for periodic timeseries when only intermediate values are known.\cr
#' For each couple of value with an increasing or deareasing segment of 
#' the sinusoid function, it is possible to estimate a minimum and maximum 
#' values using analytical algebra.\cr
#' Then the average and standard deviations of all minima and maxima are evaluated.\cr
#' It should be noted that any extremum can be estimated at least twice, one by
#' incresing segment and one by decreasing segment.
#' @examples
#' \dontrun{
#' # Generate a timeserie of time
#' time.obs <- NULL
#' for (i in 0:9) time.obs <- c(time.obs, c(0, 6, 12, 18)+i*24)
#' # For these time, generate a timeseries of temperatures
#' temp.obs <- rep(NA, length(time.obs))
#' temp.obs[3+(0:9)*4] <- rnorm(10, 25, 3)
#' temp.obs[1+(0:9)*4] <- rnorm(10, 10, 3)
#' for (i in 1:(length(time.obs)-1)) 
#'   if (is.na(temp.obs[i])) 
#'   temp.obs[i] <- mean(c(temp.obs[i-1], temp.obs[i+1]))
#'   if (is.na(temp.obs[length(time.obs)])) 
#'   temp.obs[length(time.obs)] <- temp.obs[length(time.obs)-1]/2
#' 
#' # Search for the minimum and maximum values
#' r <- minmax.periodic(time.minmax.daily=c(Min=2, Max=15), 
#' time.obs=time.obs, temp.obs=temp.obs, period=24)
#' 
#' # Estimate all the temperatures for these values
#' t <- temperature.periodic(minmax=r)
#' 
#' plot_errbar(x=t[,"time"], y=t[,"temperature"],
#' errbar.y=ifelse(is.na(t[,"sd"]), 0, 2*t[,"sd"]),
#' type="l", las=1, bty="n", errbar.y.polygon = TRUE, 
#' xlab="hours", ylab="Temperatures", ylim=c(0, 35), 
#' errbar.y.polygon.list = list(col="grey"))
#' 
#' plot_add(x=t[,"time"], y=t[,"temperature"], type="l")
#' }
#' @export

minmax.periodic <- function(time.minmax.daily, time.obs, temp.obs, period=24) {
  
maxday <- max(time.obs)/period+2
time.minmax <- NULL
for (i in 1:maxday) time.minmax <- c(time.minmax, period*(i-1)+time.minmax.daily["Min"], period*(i-1)+time.minmax.daily["Max"])

time.minmax <- time.minmax[1:which(time.minmax>max(time.obs))[1]]

res <- NULL
for (i in 1:(length(time.minmax)-1)) {
#  print(paste(1, i,length(res)))
  lim <- (time.obs>=time.minmax[i] & time.obs<=time.minmax[i+1])
  wlim <- which(lim)
  if (length(wlim)>=2) {
    couple <- t(combn(wlim, 2))
    for (j in 1:nrow(couple)) {
    x0 <- time.minmax[i]
    x1 <- time.minmax[i+1]
    x <- time.obs[couple[j, 1]]
    y <- temp.obs[couple[j, 1]]
    xp <- time.obs[couple[j, 2]]
    yp <- temp.obs[couple[j, 2]]
    C <- cos(x*pi/(x1-x0)-x0*pi/(x1-x0))
    Cp <- cos(xp*pi/(x1-x0)-x0*pi/(x1-x0))
    L <- 1/2*(1-C)
    K <- 1/2*(Cp+1) 
    y0 <- (yp/K -y/(2*K*L) +(Cp*y)/(2*K*L)) / 
      (1 -C/(4*K*L) -1/(4*K*L) +(Cp*C)/(4*K*L) +(Cp)/(4*K*L))
    y1 <- y/L -C*y0/(2*L) -y0/(2*L)
    
    if (all(time.obs!=x0)) res <- c(res, list(c(x0=unname(x0), y0=unname(y0))))
    if (all(time.obs!=x1)) res <- c(res, list(c(x1=unname(x1), y1=unname(y1))))
    }
  }
#    print(paste(2, i,length(res)))
}

names(res) <- paste0("X", 1:length(res))

res <- as.data.frame(t(as.data.frame(res)))

mxmn <- aggregate(y0 ~ x0, res, mean)
mxmn2 <- aggregate(y0 ~ x0, res, sd)

mxmn <- cbind(mxmn, mxmn2[,2])

names(mxmn) <- c("time", "temperature", "sd")

return(mxmn)

}


