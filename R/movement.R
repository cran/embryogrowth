#' movement is a function that permits to analyze movement datalogger
#' @title Analyze movement recorded within a nest with an accelerometer datalogger
#' @author Marc Girondot
#' @return The function will return a data.frame
#' @param x A data.frame with 4 columns, one for time and three for x, y, and z position
#' @param col.time Name of the column with time
#' @param col.x Name of the column with x positions
#' @param col.y Name of the column with y positions
#' @param col.z Name of the column with z positions
#' @param DaysQuiet Number of days in quiet period
#' @param SkipDays Number of days to skip before being in quiet mode
#' @param k Factor to multiply SD to prevent false positive detection
#' @param Windowsize Number of records used for moving average
#' @description This function is used to evaluate significant movement within a nest.\cr
#' The "quiet" period is the period without any expected move. It is used as a reference 
#' to detect the period with significant movements.\cr
#' It returns a data.frame with the columns:\cr
#' "Time", "x", "y", "z", \cr
#' "mvt", "mvt_standardized", "peakmvt", "running", "mvt_MA_standardized", \cr
#' "mvt_2", "mvt_2_standardized", "peakmvt_2", "running_2", "mvt_2_MA_standardized".\cr
#' mvt and mvt_2 are two different methods. Often mvt_2 is better to descrobe movments.\cr
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


movement <- function(x=stop("data.frame must be provided")  , 
                     col.time="Time"                        , 
                     col.x="x"                              , 
                     col.y="y"                              , 
                     col.z="z"                              , 
                     DaysQuiet = 40                         , 
                     SkipDays = 1                           ,
                     k=4                                    , 
                     Windowsize=15                          ) {
  data_mvtsb <- x
  
  xp <- data_mvtsb[, col.x]
  yp <- data_mvtsb[, col.y]
  zp <- data_mvtsb[, col.z]
  
  NumberRecordBeforeEmergence <- min(which((data_mvtsb[, col.time] - data_mvtsb[1, col.time])/86400 > DaysQuiet))
  if (SkipDays != 0) {
    NumberRecordsInitial <- max(which((data_mvtsb[, col.time] - data_mvtsb[1, col.time])/86400 < SkipDays))
  } else {
    NumberRecordsInitial <- 1
  }
  
  kx <- xp[-1] - xp[-length(xp)]
  kx <- c(kx[1], kx)
  ky <- yp[-1] - yp[-length(yp)]
  ky <- c(ky[1], ky)
  kz <- zp[-1] - zp[-length(zp)]
  kz <- c(kz[1], kz)
  
  mvt <- abs(kx)+abs(ky)+abs(kz)
  mvt_st <- (mvt-min(mvt, na.rm = TRUE))/(max(mvt, na.rm = TRUE)-min(mvt, na.rm = TRUE))
  mvt_st <- MovingWindow(mvt_st, window = Windowsize, hole="bothL", FUN=mean)
  
  mvt[1:(NumberRecordBeforeEmergence-1)] <- NA
  data_mvtsb <- cbind(data_mvtsb, mvt=mvt)
  data_mvtsb <- cbind(data_mvtsb, mvt_standardized=(mvt-min(mvt, na.rm = TRUE))/(max(mvt, na.rm = TRUE)-min(mvt, na.rm = TRUE)))
  data_mvtsb <- cbind(data_mvtsb, 
                      mvt_MA_standardized=MovingWindow(data_mvtsb[, "mvt_standardized"], window = Windowsize, hole="bothL", FUN=mean)
  )
  
  nomvt_mean <- mean(mvt_st[NumberRecordsInitial:(NumberRecordBeforeEmergence-1)], na.rm=TRUE)
  nomvt_sd <- sd(mvt_st[NumberRecordsInitial:(NumberRecordBeforeEmergence-1)], na.rm=TRUE)
  
  if (any(mvt_st[NumberRecordsInitial:(NumberRecordBeforeEmergence-1)] > nomvt_mean + k * nomvt_sd)) {
    inter_mvt <- mvt_st[NumberRecordsInitial:(NumberRecordBeforeEmergence-1)]
    inter_mvt[mvt_st[NumberRecordsInitial:(NumberRecordBeforeEmergence-1)] > nomvt_mean + k * nomvt_sd] <- NA
    nomvt_mean <- mean(inter_mvt, na.rm=TRUE)
    nomvt_sd <- sd(inter_mvt, na.rm=TRUE)
  }
  
  
  data_mvtsb <- cbind(data_mvtsb, peakmvt=NA)
  nvtp_peak <- data_mvtsb[, "mvt_MA_standardized"]
  data_mvtsb[which.max(nvtp_peak), "peakmvt"] <- 1
  
  diffpo <- c(0, (nvtp_peak[2:length(nvtp_peak)] - nvtp_peak[1:(length(nvtp_peak)-1)]))
  diffpo_neg <- which(diffpo < 0)
  diffpo_pos <- which(diffpo > 0)
  
  nvtp_peak[max(diffpo_neg[diffpo_neg < which.max(nvtp_peak)]):min(diffpo_pos[diffpo_pos > which.max(nvtp_peak)])] <- -1
  data_mvtsb[which.max(nvtp_peak), "peakmvt"] <- 2
  
  data_mvtsb <- cbind(data_mvtsb, running=NA)
  data_mvtsb[which(data_mvtsb[, "mvt_MA_standardized"] > nomvt_mean + k * nomvt_sd), "running"] <- 1
  
  
  # Ancienne méthode
  
  kx <- abs(xp - mean(xp[NumberRecordsInitial:NumberRecordBeforeEmergence]))
  kx <- MovingWindow(kx, window = Windowsize, hole="bothL", FUN=mean)
  kx <- kx[-1] - kx[-length(kx)]
  kx <- c(kx[1], kx)
  
  ky <- abs(yp - mean(yp[NumberRecordsInitial:NumberRecordBeforeEmergence]))
  ky <- MovingWindow(ky, window = Windowsize, hole="bothL", FUN=mean)
  ky <- ky[-1] - ky[-length(ky)]
  ky <- c(ky[1], ky)
  
  kz <- abs(zp - mean(zp[NumberRecordsInitial:NumberRecordBeforeEmergence]))
  kz <- MovingWindow(kz, window = Windowsize, hole="bothL", FUN=mean)
  kz <- kz[-1] - kz[-length(kz)]
  kz <- c(kz[1], kz)
  
  mvt <- abs(kx)+abs(ky)+abs(kz)
  mvt_st <- (mvt-min(mvt, na.rm = TRUE))/(max(mvt, na.rm = TRUE)-min(mvt, na.rm = TRUE))
  mvt_st <- MovingWindow(mvt_st, window = Windowsize, hole="bothL", FUN=mean)
  
  mvt[1:(NumberRecordBeforeEmergence-1)] <- NA
  data_mvtsb <- cbind(data_mvtsb, mvt_2=mvt)
  data_mvtsb <- cbind(data_mvtsb, mvt_2_standardized=(mvt-min(mvt, na.rm = TRUE))/(max(mvt, na.rm = TRUE)-min(mvt, na.rm = TRUE)))
  data_mvtsb <- cbind(data_mvtsb, 
                      mvt_2_MA_standardized=MovingWindow(data_mvtsb[, "mvt_2_standardized"], window = Windowsize, hole="bothL", FUN=mean)
  )
  nomvt_mean <- mean(mvt_st[NumberRecordsInitial:(NumberRecordBeforeEmergence-1)], na.rm=TRUE)
  nomvt_sd <- sd(mvt_st[NumberRecordsInitial:(NumberRecordBeforeEmergence-1)], na.rm=TRUE)
  
  if (any(mvt_st[NumberRecordsInitial:(NumberRecordBeforeEmergence-1)] > nomvt_mean + k * nomvt_sd)) {
    inter_mvt <- mvt_st[NumberRecordsInitial:(NumberRecordBeforeEmergence-1)]
    inter_mvt[mvt_st[NumberRecordsInitial:(NumberRecordBeforeEmergence-1)] > nomvt_mean + k * nomvt_sd] <- NA
    nomvt_mean <- mean(inter_mvt, na.rm=TRUE)
    nomvt_sd <- sd(inter_mvt, na.rm=TRUE)
  }
  
  data_mvtsb <- cbind(data_mvtsb, peakmvt_2=NA)
  nvtp_peak <- data_mvtsb[, "mvt_2_MA_standardized"]
  data_mvtsb[which.max(nvtp_peak), "peakmvt_2"] <- 1
  diffpo <- c(0, (nvtp_peak[2:length(nvtp_peak)] - nvtp_peak[1:(length(nvtp_peak)-1)]))
  diffpo_neg <- which(diffpo < 0)
  diffpo_pos <- which(diffpo > 0)
  
  nvtp_peak[max(diffpo_neg[diffpo_neg < which.max(nvtp_peak)]):min(diffpo_pos[diffpo_pos > which.max(nvtp_peak)])] <- -1
  data_mvtsb[which.max(nvtp_peak), "peakmvt_2"] <- 2
  
  data_mvtsb <- cbind(data_mvtsb, running_2=NA)
  data_mvtsb[which(data_mvtsb[, "mvt_2_MA_standardized"] > nomvt_mean + k * nomvt_sd), "running_2"] <- 1
  
  
  
  return(data_mvtsb)
  
}

