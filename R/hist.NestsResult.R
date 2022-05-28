#' hist.NestsResult shows the histogram of temperatures with set of nests
#' @title Show the histogram of temperatures with set of nests
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @return A list with an histogram object with information on histogram or 
#' NULL if no series was selected and the complete set of temperatures used.
#' @param x Results obtained after searchR
#' @param series Series to be used, logical (TRUE ou FALSE), numbers or names. If "all", all series are used.
#' @param ... Parameters used by hist function (example main="Title")
#' @description Show the histogram of temperatures with set of nests
#' hist(data)
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(resultNest_4p_SSM)
#' h <- hist(resultNest_4p_SSM, series=c(1:5))
#' }
#' @method hist NestsResult
#' @export

hist.NestsResult <- function(x, series="all", ...) {

# Je prends seulement les donnees de temperature et j'envoie a hist.Nests
# Ce sera plus simple pour faire les mises a jour - 30/7/2012

p3p <- list(...)
  
# j'ai un objet de resultat
# je prends les donnees
nids <- x$data
nids <- addS3Class(nids, "Nests")

L <- modifyList(list(x=nids), p3p)
L <- modifyList(L, list(series=series))

a <- do.call(hist, L)

return(invisible(a))

}
