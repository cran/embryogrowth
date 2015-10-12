#' Database of thermosensitive period of development for sex determination
#' @title Database of thermosensitive period of development for sex determination
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name TSP.list
#' @description Database of thermosensitive period of development for sex determination
#' @references Pieau, C., Dorizzi, M., 1981. Determination of temperature sensitive 
#' stages for sexual differentiation of the gonads in embryos of the turtle, 
#' Emys orbicularis. Journal of Morphology 170, 373-382.
#' Yntema, C.L., Mrosovsky, N., 1982. Critical periods and pivotal temperatures for 
#' sexual differentiation in loggerhead sea turtles. Canadian Journal of 
#' Zoology-Revue Canadienne de Zoologie 60, 1012-1016.
#' @keywords datasets
#' @family Functions for temperature-dependent sex determination
#' @usage TSP.list
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(TSP.list)
#' names(TSP.list)
#' TSP.list[["Emys_orbicularis.mass"]]
#' attributes(TSP.list[["Emys_orbicularis.mass"]])$TSP.begin.stages
#' attributes(TSP.list[["Emys_orbicularis.mass"]])$TSP.end.stages
#' }
#' @format A list with dataframes including attributes
NULL
