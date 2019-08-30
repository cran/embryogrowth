#' summary.Nests Summarize the information from a Nests object
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
#' @method summary Nests
#' @export



summary.Nests <- function(object, ...) {

	cat(paste("Number of timeseries: ", object$IndiceT["NbTS"], "\n", sep=""))
	for (i in 1:(object$IndiceT["NbTS"])) {
		nm <- paste(names(object)[i], "               " , sep="")
		pt <- paste(substr(nm, 1, 15), ":", sprintf("%.3f", max(object[[i]][,1])/1440), " days\n", sep="")
		cat(pt)
	}

}
