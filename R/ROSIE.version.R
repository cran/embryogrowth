#' Database of information for incubation of turtles
#' @title Version of database of TSD information for reptiles
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @return The date of the latest updated version
#' @description Return the date of the most recent update of the ROSIE database.
#' @family Functions for temperature-dependent sex determination
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' ROSIE.version()
#' }
#' @export

ROSIE.version <- function() {
  # data("DatabaseTSD", package = "embryogrowth")
  cat("This is a modified version of the ROSIE database.")
  return("1.0.3 modified 2025-06-03")
}
