#' web.tsd runs a shiny application for basic functions of tsd function
#' @title Run a shiny application for basic functions of tsd function
#' @author Marc Girondot
#' @return Nothing
#' @description Run a shiny application for basic functions of tsd function.
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' web.tsd()
#' }
#' @export


web.tsd <- function() {
  
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is absent; Please install it first")
  }
  
getFromNamespace("runApp", ns="shiny")(appDir = system.file("shiny", package="embryogrowth"), 
                                       launch.browser =TRUE)

}
