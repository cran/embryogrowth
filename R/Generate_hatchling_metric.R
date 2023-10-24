#' Generate_hatchling_metric Generate a data.frame that can be used as hatchling.metric value for searchR()
#' @title Generate a data.frame that can be used as hatchling.metric value for searchR()
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A data.frame with size or mass at hatching for each nest
#' @param series Name of series or object from searchR()
#' @param hatchling.metric Size or mass at hatching. Will be recycled if necessary
#' @param previous Previous formated hatchling.metric data
#' @description Generate a data.frame that can be used as hatchling.metric value for searchR()
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(resultNest_4p_SSM)
#' testsize1 <- Generate_hatchling_metric(resultNest_4p_SSM)
#' testsize2 <- Generate_hatchling_metric(series=resultNest_4p_SSM,  
#' 	hatchling.metric=c(Mean=39.3, SD=1.92))
#' }
#' @export


Generate_hatchling_metric <-
function(series=stop("A result object or names of series must be provided"), 
         hatchling.metric=NULL, previous=NULL) {

	if (is.null(hatchling.metric) & (!inherits(series, "NestsResult"))) { #(class(series)!="NestsResult")) {
		stop("hatchling.metric or a result from searchR() must be provided")
	}
	
	if (inherits(series, "NestsResult")) { #(class(series)=="NestsResult") {
		if (is.null(hatchling.metric)) testec <- series$hatchling.metric
		series <- series$data
	}
	
	if (inherits(series, "Nests")) { #(class(series)=="Nests") {
		series <- names(series)
		series <- series[1:(series$IndiceT["NbTS"])]
	}
	
	
	if (!is.null(hatchling.metric)) {
		if (is.na(hatchling.metric["Mean"]) | is.na(hatchling.metric["SD"]) | length(hatchling.metric["Mean"])!=length(hatchling.metric["Mean"])) {
			return("hatchling.metric must be a vector with same number of Mean and SD values")
		} else {
			mean <- rep(hatchling.metric["Mean"], length(series))[1:length(series)]
			sd <- rep(hatchling.metric["SD"], length(series))[1:length(series)]
			hatchling.metric_ec <- data.frame(Mean=mean, SD=sd, row.names=series)
		}
	}
	
	if (!is.null(previous)) {
	  hatchling.metric <- rbind(previous, hatchling.metric_ec)
	} else {
	  hatchling.metric <- hatchling.metric_ec
	}
		
	return(hatchling.metric)	

}
