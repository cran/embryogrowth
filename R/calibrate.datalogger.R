#' calibrate.datalogger calibrates data loggers and correct time series of temperatures.
#' @title Calibrate data loggers and correct time series of temperatures
#' @author Marc Girondot
#' @return The function will return a corrected time series of temperatures as a vector if se.fit is FALSE or a list if se.fit is TRUE.
#' @param control.temperatures The true temperatures during the calibration process
#' @param read.temperatures The read temperatures during the calibration process
#' @param temperatures.series The temperatures to be converted using calibration
#' @param gam Does gam should be used (TRUE) or glm (FALSE).
#' @param se.fit Do standard errors are to be returned
#' @description Calibrate a time series of temperatures. Use or gam or glm. If
#' no temperatures.series is given, it will use the read.temperatures.
#' @references
#' \insertRef{11124}{embryogrowth}\cr
#' @family Data loggers utilities
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' calibrate.datalogger(control.temperatures=20:30, 
#'                      read.temperatures=(20:30)+rnorm(11))
#' }
#' @export


calibrate.datalogger <- function(
                     control.temperatures=stop("Control temperatures is missing"), 
                     read.temperatures=stop("Read temperatures must be indicated"), 
                     temperatures.series=NULL, 
                     gam=TRUE, se.fit=TRUE) {
  
  if ((!requireNamespace("gam", quietly = TRUE)) | !gam) {
    g <- glm(control.temperatures ~ read.temperatures)
  } else {
    g <- getFromNamespace("gam", ns="gam")(control.temperatures ~ read.temperatures)
  }
  if (is.null(temperatures.series)) temperatures.series <- read.temperatures
    return(predict(g, data.frame(read.temperatures=temperatures.series), se.fit=se.fit))
}

