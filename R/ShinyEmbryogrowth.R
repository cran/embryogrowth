#' ShinyEmbryogrowth runs a shiny application for basic functions of embryogrowth
#' @title Run a shiny application for basic functions of embryogrowth
#' @author Marc Girondot
#' @author Maria Sousa Martins
#' @return Nothing
#' @description Run a shiny application for basic functions of embryogrowth
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' ShinyEmbryogrowth()
#' }
#' @export


.ShinyEmbryogrowth <- function() {

runApp(appDir = system.file("shiny", package="embryogrowth"),launch.browser =TRUE)

}
