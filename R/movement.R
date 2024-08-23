#' movement is a function that permits to analyze movement datalogger
#' @title Analyze movement recorded within a nest with an accelerometer datalogger
#' @author Marc Girondot
#' @return The function will return a list
#' @param x A data.frame with 4 columns, one for time and three for x, y, and z position
#' @param col.time Name of the column with time
#' @param col.x Name of the column with x positions
#' @param col.y Name of the column with y positions
#' @param col.z Name of the column with z positions
#' @param NumberRecordBeforeEmergence Number of records in quiet period
#' @param k Factor to multiply SD to prevent false positive detection
#' @param Windowsize Number of records used for moving average
#' @description This function is used to evaluate significant movement within a nest.
#' @references
#' \insertRef{11893}{embryogrowth}\cr
#' @family Data loggers utilities
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' mv <- movement(x=dataf, 
#'                col.time="Time", 
#'                col.x="x", col.y="y", col.z="z")
#' }
#' @export


movement <- function(x=stop("data.frame must be provided"), 
                     col.time="Time", 
                     col.x="x", col.y="y", col.z="z", 
                     NumberRecordBeforeEmergence = 1900, 
                     k=4, Windowsize=15) {
  data_mvtsb <- x
  
  xp <- data_mvtsb[, col.x]
  yp <- data_mvtsb[, col.y]
  zp <- data_mvtsb[, col.z]
  
  kx <- abs(xp - mean(xp[1:NumberRecordBeforeEmergence]))
  kx <- MovingWindow(kx, window = Windowsize, hole="bothL", FUN=mean)
  kx <- kx[-1] - kx[-length(kx)]
  kx <- c(kx[1], kx)
  
  ky <- abs(yp - mean(yp[1:NumberRecordBeforeEmergence]))
  ky <- MovingWindow(ky, window = Windowsize, hole="bothL", FUN=mean)
  ky <- ky[-1] - ky[-length(ky)]
  ky <- c(ky[1], ky)
  
  kz <- abs(zp - mean(zp[1:NumberRecordBeforeEmergence]))
  kz <- MovingWindow(kz, window = Windowsize, hole="bothL", FUN=mean)
  kz <- kz[-1] - kz[-length(kz)]
  kz <- c(kz[1], kz)
  
  mvt <- abs(kx)+abs(ky)+abs(kz)
  
  data_mvtsb <- cbind(data_mvtsb, mvt=(mvt-min(mvt))/(max(mvt)-min(mvt)))
  
  peakmvt <- data_mvtsb[which.max(data_mvtsb[, "mvt"]), ]
  
  mvtrunning <- data_mvtsb[data_mvtsb[, "mvt"] > mean(data_mvtsb$mvt[1:NumberRecordBeforeEmergence])+
                             k*sd(data_mvtsb$mvt[1:NumberRecordBeforeEmergence]), ]
  
  return(list(index.peak=which.max(mvtrunning$mvt), 
         time.peak=peakmvt[, col.time], 
         movements=mvtrunning))
  
}

