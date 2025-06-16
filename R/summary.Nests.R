#' summary.Nests2 Summarize the information from a Nests object
#' @title Summarize the information from a Nests object.
#' @author Marc Girondot
#' @return None
#' @param object A object obtained after FormatNests()
#' @param ... Not used
#' @description Summarize the information from a Nests object:\cr
#' - Name of the nests, total incubation length and average temperature
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest, previous=NULL)
#' summary(formated)
#' }
#' @method summary Nests2
#' @export



summary.Nests2 <- function(object, ...) {
  
  cat(paste("Number of timeseries: ", object$IndiceT["NbTS"], "\n", sep=""))
  for (i in object$Names) {
    cat(paste0(substr(paste0(i, "               "), 1, 15), ":       ", sprintf("%.3f", max(object$Nests[[i]]$data[,1])/1440), " days\n", 
                 ifelse(!is.null(object$Nests[[i]]$Informations), paste0("Informations: ", object$Nests[[i]]$Informations, "\n"), ""), 
                 ifelse(!is.null(object$Nests[[i]]$Longitude), paste0("Longitude:             ", object$Nests[[i]]$Longitude, "\n"), ""), 
                 ifelse(!is.null(object$Nests[[i]]$Latitude), paste0("Latitude:              ", object$Nests[[i]]$Latitude, "\n"), ""), 
                 ifelse(!is.null(object$Nests[[i]]$hatchling.metric.mean), paste0("Hatchling metric mean: ", object$Nests[[i]]$hatchling.metric.mean, "\n"), ""), 
                 ifelse(!is.null(object$Nests[[i]]$hatchling.metric.sd), paste0("Hatchling metric sd:   ", object$Nests[[i]]$hatchling.metric.sd, "\n"), ""), 
                 ifelse(!is.null(object$Nests[[i]]$Males), paste0("Number of males:                 ", object$Nests[[i]]$Males, "\n"), ""), 
                 ifelse(!is.null(object$Nests[[i]]$Females), paste0("Number of females:                 ", object$Nests[[i]]$Females, "\n"), "")))
    cat("________________________________\n")
  }
  
}
