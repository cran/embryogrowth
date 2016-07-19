#' plotR_hist shows the histogram of temperatures with set of nests and the R function superimpose
#' @title Shows the histogram of temperatures with set of nests and the R function superimpose
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return A list with data.frame with the confidence interval and the average.
#' @param result Result data
#' @param ... Parameters used by hist or plotR functions
#' @param ylimH Scale of histogram using ylimH=c(min, max)
#' @param atH Position of ticks for scale of histogram
#' @param ylabH Label for histogram scale
#' @description Shows the histogram of temperatures with set of nests and the R function superimpose
#' plotR_hist(data)
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(resultNest_4p)
#' plotR_hist(resultNest_4p)
#' plotR_hist(resultNest_4p, ylim=c(0,0.3), ylimH=c(0,0.5), atH=c(0, 0.1, 0.2))
#' }
#' @export

plotR_hist <- function(result=NULL, ..., ylimH=NULL, atH=NULL, ylabH="Frequency of temperatures") {

def.par <- par(no.readonly = TRUE) # save default, for resetting...

p3p <- list(...)
  
  
if (class(result)=="list") {
	stop("Only one series can be used with plotR_hist()")
}

nids <- result

par(mar = c(def.par[["mar"]][1:3], 5.1))

L <- modifyList(list(result=result), p3p)
lp <- list("series")
L2 <- L[!(names(L) %in% lp)]


output <- do.call(plotR, L2) 

ylim <- ScalePreviousPlot()$ylim[c("begin", "end")]

par(new=TRUE)

L <- modifyList(list(x=result), p3p)

L["ylim"] <- NULL
L["SE"] <- NULL
L["legend"] <- NULL

if (!is.null(ylimH)) L <- c(L, list(ylim=ylimH))

lp <- list("show.box", "set.par", "parameters", "fixed.parameters", "legend", "size", 
"scaleY", "replicate.CI", "ylabH", "xlimR", "ltyCI", "lwdCI", "xlimSE", "x.SE")
L <- L[!(names(L) %in% lp)]

# x2 <- (par("usr")[1]+par("usr")[2]*26)/27
# x1 <- x2*26-par("usr")[2]/0.04
xlim <- ScalePreviousPlot()$xlim[c("begin", "end")]

L <- modifyList(L, list(xlab="", ylab="", main="", axes=FALSE, freq=FALSE, 
	xlim=xlim)) 

a <- do.call(hist, L) 


axis(side=4, ylim=par("yaxp")[1:2], las=1, at=atH, 
     cex.axis=ifelse(any(names(p3p)=="cex.axis"), p3p$cex.axis, par("cex")))
mtext(ylabH, side=4, line=3, 
      cex=ifelse(any(names(p3p)=="cex.lab"), p3p$cex.lab, par("cex")))
par(new=TRUE)
# je retablis l'echelle des y et celle de R
plot(x = 1, y=1, ylim=ylim, xlim=xlim, xlab="", ylab="", axes=FALSE, bty="n", type="n")

# par(def.par)  #- reset to default

return(invisible(output))
}
