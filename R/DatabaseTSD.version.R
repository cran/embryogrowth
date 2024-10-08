#' Database of information for incubation of turtles
#' @title Version of database of TSD information for reptiles
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @return The date of the lastest updated version
#' @description Return the date of the most recent update of the database.
#' @family Functions for temperature-dependent sex determination
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' DatabaseTSD.version()
#' }
#' @export

DatabaseTSD.version <- function() {
  # data("DatabaseTSD", package = "embryogrowth")
  return(max(embryogrowth::DatabaseTSD$Version)[1])
}
