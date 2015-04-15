#' plotR_hist shows the histogram of temperatures with set of nests and the R function superimpose
#' @title Shows the histogram of temperatures with set of nests and the R function superimpose
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return Nothing
#' @param x Result data
#' @param ... Parameters used by hist or plotR functions
#' @param ylimH Scale of histogram using ylimH=c(min, max)
#' @param ylabH Label for histogram scale
#' @description Shows the histogram of temperatures with set of nests and the R function superimpose
#' plotR_hist(data)
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(resultNest_4p)
#' plotR_hist(resultNest_4p)
#' }
#' @export

plotR_hist <- function(x, ..., ylimH=NULL, ylabH="Frequency of temperatures") {

def.par <- par(no.readonly = TRUE) # save default, for resetting...
  
  
if (class(x)=="list") {

	print("Only one series can be used with plotR_hist()")

} else {

nids <- x

par(mar = c(def.par[["mar"]][1:3], 5.1))

L <- list(result=x, ...)

a <- do.call(plotR, L) 

par(new=TRUE)

L <- list(x=x, ...)

L["ylim"] <- NULL
L["SE"] <- NULL
L["legend"] <- NULL

if (!is.null(ylimH)) L <- c(L, list(ylim=ylimH))

lp <- list("show.box", "set.par", "parameters", "fixed.parameters", "legend", "size", 
"scaleY", "replicate.CI", "ylabH", "xlimR", "ltyCI", "lwdCI")
L <- L[!(names(L) %in% lp)]


#L <- ifelse(length(which(names(L)=="show.box"))==0, 0, L[-which(names(L)=="show.box")])
#L <- ifelse(length(which(names(L)=="set.par"))==0, 0, L[-which(names(L)=="set.par")])
#L <- ifelse(length(which(names(L)=="parameters"))==0, 0, L[-which(names(L)=="parameters")])
#L <- ifelse(length(which(names(L)=="fixed.parameters"))==0, L, L[-which(names(L)=="fixed.parameters")])
#L <- ifelse(length(which(names(L)=="legend"))==0, L, L[-which(names(L)=="legend")])
#L <- ifelse(length(which(names(L)=="size"))==0, L, L[-which(names(L)=="size")])
#L <- ifelse(length(which(names(L)=="scaleY"))==0, L, L[-which(names(L)=="scaleY")])
#L <- ifelse(length(which(names(L)=="replicate.CI"))==0, L, L[-which(names(L)=="replicate.CI")])
#L <- ifelse(length(which(names(L)=="ylabH"))==0, L, L[-which(names(L)=="ylabH")])


x2 <- (par("usr")[1]+par("usr")[2]*26)/27
x1 <- x2*26-par("usr")[2]/0.04

L <- modifyList(L, list(xlab="", ylab="", main="", axes=FALSE, freq=FALSE, 
	xlim=c(x1, x2))) 

a <- do.call(hist, L) 

axis(side=4, ylim=par("yaxp")[1:2], las=1)
mtext(ylabH, side=4, line=3)

# par(def.par)  #- reset to default


}
}
