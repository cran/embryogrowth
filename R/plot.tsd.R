#' plot.tsd plot result of tsd() that best describe temperature-dependent sex determination
#' @title Plot results of tsd() that best describe temperature-dependent sex determination
#' @author Marc Girondot
#' @return Nothing
#' @param x A result file generated by tsd()
#' @param ... Parameters for plot()
#' @param show.observations Should the observations be shown
#' @param show.observations.sd Should the observations SD be shown
#' @param show.model Should the model be shown
#' @param show.PTRT Should the P and TRT information be shown
#' @param resultmcmc A result of tsd_MHmcmc()
#' @param chain What chain to be used is resultmcmc is provided
#' @param temperatures.plot Temperatures used for showing curves of sex ratio
#' @param durations.plot Durations used for showing curves of sex ratio
#' @param replicate.CI replicate.CI replicates from the hessian matrix to estimate CI 
#' @param range.CI The range of confidence interval for estimation, default=0.95
#' @param l Sex ratio limits to define TRT are l and 1-l
#' @param las.x las parameter for x axis
#' @param las.y las parameter for y axis
#' @param lab.PT Label to describe pivotal temperature
#' @param lab.TRT Label to describe transitional range of temperature
#' @param males.freq Should the graph uses males relative frequency [TRUE] or females [FALSE]
#' @param mar The par("mar") parameter
#' @param col.TRT The color of TRT
#' @param col.TRT.CI The color of CI of TRT based on range.CI
#' @param col.PT.CI The color of CI of PT based on range.CI
#' @param show.CI Do the CI for the curve should be shown
#' @param warn Do the warnings must be shown ? TRUE or FALSE
#' @param use.ggplot Use ggplot graphics (experimental). TRUE or FALSE
#' @description Plot the estimates that best describe temperature-dependent sex determination.\cr
#' \insertRef{1515}{embryogrowth}\cr
#' \insertRef{3534}{embryogrowth}\cr
#' \insertRef{11754}{embryogrowth}\cr
#' \insertRef{5790}{embryogrowth}\cr
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' CC_AtlanticSW <- subset(DatabaseTSD, RMU.2010=="Atlantic, SW" & 
#'                           Species=="Caretta caretta" & (!is.na(Sexed) & Sexed!=0))
#' tsdL <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature.set, 
#'                                  equation="logistic"))
#' # By default, it will return a ggplot object
#' # Here I show the advantage of using a ggplot object
#' g <- plot(tsdL)
#' # You can remove named layers. For example:
#' names(g$layers)
#' g$layers["Observations"] <- NULL; plot(g)
#' # And add some
#' # Due to a bug in ggplot, it is necessary to remove all names to obtain correct legends
#' names(g$layers) <- NULL
#' g + geom_point(data=CC_AtlanticSW, aes(x=Incubation.temperature.set, y=Males/Sexed, 
#'                size = Sexed), inherit.aes = FALSE, show.legend = TRUE, shape=19)
#' # Force to use the original plot
#' plot(tsdL, use.ggplot = FALSE)
#' }
#' @family Functions for temperature-dependent sex determination
#' @method plot tsd
#' @export



plot.tsd <- function(x, ..., 
                     show.observations=TRUE                                                , 
                     show.observations.sd=TRUE                                             ,
                     show.model=TRUE                                                       , 
                     males.freq=TRUE                                                       , 
                     show.PTRT=TRUE                                                        , 
                     las.x=1                                                               , 
                     las.y=1                                                               , 
                     lab.PT=paste0("Pivotal ", x$type)                                      , 
                     resultmcmc = NULL                                                     ,
                     chain=1                                                               , 
                     l=0.05                                                                , 
                     replicate.CI=10000                                                    , 
                     range.CI=0.95                                                         , 
                     mar=c(4, 4, 4, 1)+0.4                                                 ,
                     temperatures.plot=seq(from=25, to=35, by=0.1)                         , 
                     durations.plot=seq(from=40, to=70, by=0.1)                            , 
                     lab.TRT=paste0("Transitional range of ",  x$type, "s l=",x$l*100,"%") , 
                     col.TRT="gray"                                                        , 
                     col.TRT.CI=rgb(0.8, 0.8, 0.8, 0.8)                                    , 
                     col.PT.CI=rgb(0.8, 0.8, 0.8, 0.8)                                     , 
                     show.CI=TRUE                                                          , 
                     warn = TRUE                                                           , 
                     use.ggplot = TRUE                                                     ) {
  
  # show.observations=TRUE; males.freq=TRUE; show.PTRT = TRUE; las.x=1; las.y=1; lab.PT=paste("Pivotal ", x$type); resultmcmc = NULL; chain=1; l=0.05; replicate.CI=10000; temperatures.plot=seq(from=20, to=35, by=0.1); durations.plot=seq(from=40, to=70, by=0.1);TRT.limits=c(9, 90); precision=15; range.CI=0.95; mar=c(4, 4, 6, 1)+0.4; lab.TRT=paste0("Transitional range of ", x$type, "s l=",x$l*100,"%"); col.TRT="gray"; col.TRT.CI=rgb(0.8, 0.8, 0.8, 0.5); col.PT.CI=rgb(0.8, 0.8, 0.8, 0.5); show.CI=TRUE; warn=TRUE
  
  L <- list(...) # L <- list()
  
  if (use.ggplot & all(rownames(installed.packages()) != "ggplot2")) {
    warning("Install ggplot2 package to use the option use.ggplot = TRUE")
    use.ggplot = FALSE
  }
  
  males <- x$males
  females <- x$females
  N <- x$N
  temperatures <- x$temperatures
  equation <- x$equation  
  
  if (x$type != "temperature") temperatures.plot <- durations.plot
  
  xlll <- ifelse(x$type=="temperature", expression("Temperature in " * degree * "C"), 
                 "Duration in days")
  
  if (males.freq) {
    L1 <- modifyList(list(x=temperatures, y=males/N, bty="n", type="n", xlab=xlll, 
                          ylab="Male relative frequency"), L)
  } else {
    L1 <- modifyList(list(x=temperatures, y=females/N, bty="n", type="n", xlab=xlll, 
                          ylab="Female relative frequency"), L)	
  }
  L1 <- modifyList(L1, list(ylim=c(0,1), xaxt="n", las=las.y))
  
  if (is.null(L$xlim)) {
    L1 <- modifyList(L1, list(xlim=c(floor(min(temperatures)), floor(1+max(temperatures)))))
  }
  
  temperatures.plot <- seq(from=L1$xlim[1], to=L1$xlim[2], by=0.1) 
  
  o <- P_TRT(x=x, resultmcmc=resultmcmc, chain=chain, l=l, 
             replicate.CI=replicate.CI, temperatures=temperatures.plot, 
             probs = c((1-range.CI)/2, 0.5, 1-(1-range.CI)/2), warn=warn)
  
  if (!use.ggplot) {
    L2 <- L1[!(names(L1) %in% c("errbar.tick", "errbar.lwd"
                                , "errbar.lty", "errbar.col"
                                , "errbar.y.polygon"
                                , "errbar.y.polygon.list"))]
    par(mar=mar)
    a <- do.call(plot, L2)
    x2 <- (par("usr")[1]+par("usr")[2]*26)/27
    x1 <- x2*26-par("usr")[2]/0.04
    
  } else {
    x2 <- L1$xlim[2]
    x1 <- L1$xlim[1]
    g <- ggplot(data = NULL, aes(x=NULL, y=NULL)) + ylab(L1$ylab) + xlab(L1$xlab) + 
      theme_classic() + 
      theme(plot.margin=grid::unit(c(5,5,5,5), "mm")) + theme(axis.line.y = element_blank()) +
      scale_y_continuous(limits = c(0, 1.2), expand = c(0, 0), breaks = seq(from=0, to=1, by=0.1)) + 
      scale_x_continuous(limits = L1$xlim, expand = c(0, 0)) + 
      geom_segment(aes(x=x1, y=0, xend = x1, yend = 1), linewidth = 1)
      names(g$layers)[length(g$layers)] <- "Base"
  }
  
  
  cex.x <- par("cex.axis")
  if (!is.null(L$cex.axis)) cex.x <- L$cex.axis
  
  ptax <- TRUE
  if (!is.null(L[["axes"]])) if (isFALSE(L[["axes"]])) ptax <- FALSE
  if (!is.null(L[["xaxt"]])) if (L[["xaxt"]] == "n") ptax <- FALSE
  
  if (!use.ggplot) {
    if (ptax) axis(1, at=x1:x2, las=las.x, cex.axis=cex.x)
    par(xpd=FALSE)
  }
  
  # je trace la TRT centree sur P
  
  
  
  if ((equation != "GSD") & (show.PTRT)) {
    if (!use.ggplot) {
      if (all(colnames(o$P_TRT_quantiles) != "PT_low")) {
        par(xpd=FALSE)
        polygon(c(o$P_TRT_quantiles[2, "lower.limit.TRT"], o$P_TRT_quantiles[2, "lower.limit.TRT"], o$P_TRT_quantiles[2, "higher.limit.TRT"], o$P_TRT_quantiles[2, "higher.limit.TRT"]), c(0,1,1,0), border=NA, col=col.TRT)  
        # limites de la limite basse de la TRT
        polygon(c(o$P_TRT_quantiles[1, "lower.limit.TRT"], o$P_TRT_quantiles[1, "lower.limit.TRT"], o$P_TRT_quantiles[3, "lower.limit.TRT"], o$P_TRT_quantiles[3, "lower.limit.TRT"]), c(0,1,1,0), border=NA, col=col.TRT.CI)
        # limites de la limite haute de la TRT
        polygon(c(o$P_TRT_quantiles[1, "higher.limit.TRT"], o$P_TRT_quantiles[1, "higher.limit.TRT"], o$P_TRT_quantiles[3, "higher.limit.TRT"], o$P_TRT_quantiles[3, "higher.limit.TRT"]), c(0,1,1,0), border=NA, col=col.TRT.CI)
        par(xpd=TRUE)
        segments(o$P_TRT_quantiles[2, "lower.limit.TRT"], 0, o$P_TRT_quantiles[2, "lower.limit.TRT"], 1.15, lty=3)
        segments(o$P_TRT_quantiles[2, "higher.limit.TRT"], 0, o$P_TRT_quantiles[2, "higher.limit.TRT"], 1.15, lty=3)
        text(x=o$P_TRT_quantiles[2, "lower.limit.TRT"]+o$P_TRT_quantiles[2, "TRT"]/2, y=1.2, lab.TRT)
        
        if (any(colnames(o$P_TRT_quantiles) == "PT")) {
        # limites de la PT
        polygon(c(o$P_TRT_quantiles[1, "PT"], o$P_TRT_quantiles[1, "PT"], o$P_TRT_quantiles[3, "PT"], o$P_TRT_quantiles[3, "PT"]), c(0,1,1,0), border=NA, col=col.PT.CI)  
        par(xpd=TRUE)
        segments(o$P_TRT_quantiles[2, "PT"], 0, o$P_TRT_quantiles[2, "PT"], 1.05, lty=4)
        text(x=o$P_TRT_quantiles[2, "PT"], y=1.1, lab.PT)
        }

        par(xpd=FALSE)
      } else {
        par(xpd=FALSE)
        polygon(c(o$P_TRT_quantiles[2, "lower.limit.TRT_low"], o$P_TRT_quantiles[2, "lower.limit.TRT_low"], o$P_TRT_quantiles[2, "higher.limit.TRT_low"], o$P_TRT_quantiles[2, "higher.limit.TRT_low"]), c(0,1,1,0), border=NA, col=col.TRT)  
        # limites de la limite basse de la TRT
        polygon(c(o$P_TRT_quantiles[1, "lower.limit.TRT_low"], o$P_TRT_quantiles[1, "lower.limit.TRT_low"], o$P_TRT_quantiles[3, "lower.limit.TRT_low"], o$P_TRT_quantiles[3, "lower.limit.TRT_low"]), c(0,1,1,0), border=NA, col=col.TRT.CI)
        # limites de la limite haute de la TRT
        polygon(c(o$P_TRT_quantiles[1, "higher.limit.TRT_low"], o$P_TRT_quantiles[1, "higher.limit.TRT_low"], o$P_TRT_quantiles[3, "higher.limit.TRT_low"], o$P_TRT_quantiles[3, "higher.limit.TRT_low"]), c(0,1,1,0), border=NA, col=col.TRT.CI)
        # limites de la PT
        polygon(c(o$P_TRT_quantiles[1, "PT_low"], o$P_TRT_quantiles[1, "PT_low"], o$P_TRT_quantiles[3, "PT_low"], o$P_TRT_quantiles[3, "PT_low"]), c(0,1,1,0), border=NA, col=col.PT.CI)  
        par(xpd=TRUE)
        segments(o$P_TRT_quantiles[2, "PT_low"], 0, o$P_TRT_quantiles[2, "PT_low"], 1.05, lty=4)
        segments(o$P_TRT_quantiles[2, "lower.limit.TRT_low"], 0, o$P_TRT_quantiles[2, "lower.limit.TRT_low"], 1.15, lty=3)
        segments(o$P_TRT_quantiles[2, "higher.limit.TRT_low"], 0, o$P_TRT_quantiles[2, "higher.limit.TRT_low"], 1.15, lty=3)
        text(x=o$P_TRT_quantiles[2, "PT_low"], y=1.1, lab.PT)
        text(x=o$P_TRT_quantiles[2, "PT_low"], y=1.2, lab.TRT)
        par(xpd=FALSE)
        polygon(c(o$P_TRT_quantiles[2, "lower.limit.TRT_high"], o$P_TRT_quantiles[2, "lower.limit.TRT_high"], o$P_TRT_quantiles[2, "higher.limit.TRT_high"], o$P_TRT_quantiles[2, "higher.limit.TRT_high"]), c(0,1,1,0), border=NA, col=col.TRT)  
        # limites de la limite basse de la TRT
        polygon(c(o$P_TRT_quantiles[1, "lower.limit.TRT_high"], o$P_TRT_quantiles[1, "lower.limit.TRT_high"], o$P_TRT_quantiles[3, "lower.limit.TRT_high"], o$P_TRT_quantiles[3, "lower.limit.TRT_high"]), c(0,1,1,0), border=NA, col=col.TRT.CI)
        # limites de la limite haute de la TRT
        polygon(c(o$P_TRT_quantiles[1, "higher.limit.TRT_high"], o$P_TRT_quantiles[1, "higher.limit.TRT_high"], o$P_TRT_quantiles[3, "higher.limit.TRT_high"], o$P_TRT_quantiles[3, "higher.limit.TRT_high"]), c(0,1,1,0), border=NA, col=col.TRT.CI)
        # limites de la PT
        polygon(c(o$P_TRT_quantiles[1, "PT_high"], o$P_TRT_quantiles[1, "PT_high"], o$P_TRT_quantiles[3, "PT_high"], o$P_TRT_quantiles[3, "PT_high"]), c(0,1,1,0), border=NA, col=col.PT.CI)  
        par(xpd=TRUE)
        segments(o$P_TRT_quantiles[2, "PT_high"], 0, o$P_TRT_quantiles[2, "PT_high"], 1.05, lty=4)
        segments(o$P_TRT_quantiles[2, "lower.limit.TRT_high"], 0, o$P_TRT_quantiles[2, "lower.limit.TRT_high"], 1.15, lty=3)
        segments(o$P_TRT_quantiles[2, "higher.limit.TRT_high"], 0, o$P_TRT_quantiles[2, "higher.limit.TRT_high"], 1.15, lty=3)
        text(x=o$P_TRT_quantiles[2, "PT_high"], y=1.1, lab.PT)
        text(x=o$P_TRT_quantiles[2, "PT_high"], y=1.2, lab.TRT)
        par(xpd=FALSE)
      }
    } else {
      ## Je suis en ggplot
      if (all(colnames(o$P_TRT_quantiles)!="PT_low")) {
        # Patron TSD I
        g <- g + geom_rect(aes(xmin = max(x1, o$P_TRT_quantiles["50%", "lower.limit.TRT"]), 
                               ymin = 0, 
                               xmax= min(x2, o$P_TRT_quantiles["50%", "higher.limit.TRT"]), 
                               ymax = 1), colour = NA, fill = col.TRT)
        names(g$layers)[length(g$layers)] <- "TRT.median"
        # limites de la limite basse de la TRT
        g <- g + geom_rect(aes(xmin = max(x1, o$P_TRT_quantiles["2.5%", "lower.limit.TRT"]), 
                               ymin = 0, 
                               xmax= min(x2, o$P_TRT_quantiles["97.5%", "lower.limit.TRT"]), 
                               ymax = 1), colour = NA, fill = col.TRT.CI)
        names(g$layers)[length(g$layers)] <- "lower.limit.TRT.credible.interval"
        # limites de la limite haute de la TRT
        g <- g + geom_rect(aes(xmin = max(x1, o$P_TRT_quantiles["2.5%", "higher.limit.TRT"]), 
                               ymin = 0, 
                               xmax= min(x2, o$P_TRT_quantiles["97.5%", "higher.limit.TRT"]), 
                               ymax = 1), colour = NA, fill = col.TRT.CI)
        names(g$layers)[length(g$layers)] <- "higher.limit.TRT.credible.interval"
        g <- g + geom_segment(aes(x = o$P_TRT_quantiles["50%", "lower.limit.TRT"], y=0, 
                                  xend = o$P_TRT_quantiles["50%", "lower.limit.TRT"], yend=1.15), linetype = 3)
        names(g$layers)[length(g$layers)] <- "lower.limit.TRT.median"
        g <- g + geom_segment(aes(x = o$P_TRT_quantiles["50%", "higher.limit.TRT"], y=0, 
                                  xend = o$P_TRT_quantiles["50%", "higher.limit.TRT"], yend=1.15), linetype = 3)
        names(g$layers)[length(g$layers)] <- "higher.limit.TRT.median"
        
        g <- g + geom_text(aes(x=o$P_TRT_quantiles[2, "lower.limit.TRT"]+o$P_TRT_quantiles[2, "TRT"]/2, y=1.15), label = lab.TRT)
        names(g$layers)[length(g$layers)] <- "TRT.label"
        
        if (any(colnames(o$P_TRT_quantiles) == "PT")) {
        # limites de la PT
        g <- g + geom_rect(aes(xmin = o$P_TRT_quantiles["2.5%", "PT"], 
                               ymin = 0, 
                               xmax= o$P_TRT_quantiles["97.5%", "PT"], 
                               ymax = 1), colour = NA, fill = col.PT.CI)
        names(g$layers)[length(g$layers)] <- "PT.credible.interval"
        g <- g + geom_segment(aes(x = o$P_TRT_quantiles["50%", "PT"], y=0, 
                                  xend = o$P_TRT_quantiles["50%", "PT"], yend=1.05), linetype = 4)
        names(g$layers)[length(g$layers)] <- "PT.median"
        g <- g + geom_text(aes(x=o$P_TRT_quantiles["50%", "PT"], y=1.08), label = lab.PT)
        names(g$layers)[length(g$layers)] <- "PT.label"
        }
        

      } else {
        # Patron TSD II
        
        g <- g + geom_rect(aes(xmin = o$P_TRT_quantiles["50%", "lower.limit.TRT_low"], 
                               ymin = 0, 
                               xmax= o$P_TRT_quantiles["50%", "higher.limit.TRT_low"], 
                               ymax = 1), colour = NA, fill = col.TRT)
        names(g$layers)[length(g$layers)] <- "low.TRT.median"
        # limites de la limite basse de la TRT
        g <- g + geom_rect(aes(xmin = o$P_TRT_quantiles["2.5%", "lower.limit.TRT_low"], 
                               ymin = 0, 
                               xmax= o$P_TRT_quantiles["97.5%", "lower.limit.TRT_low"], 
                               ymax = 1), colour = NA, fill = col.TRT.CI)
        names(g$layers)[length(g$layers)] <- "low.TRT.lower.limit.credible.interval"
        # limites de la limite haute de la TRT
        g <- g + geom_rect(aes(xmin = o$P_TRT_quantiles["2.5%", "higher.limit.TRT_low"], 
                               ymin = 0, 
                               xmax= o$P_TRT_quantiles["97.5%", "higher.limit.TRT_low"], 
                               ymax = 1), colour = NA, fill = col.TRT.CI)
        names(g$layers)[length(g$layers)] <- "low.TRT.upper.limit.credible.interval"
        # limites de la PT
        g <- g + geom_rect(aes(xmin = o$P_TRT_quantiles["2.5%", "PT_low"], 
                               ymin = 0, 
                               xmax= o$P_TRT_quantiles["97.5%", "PT_low"], 
                               ymax = 1), colour = NA, fill = col.PT.CI)
        names(g$layers)[length(g$layers)] <- "low.PT.credible.interval"
        g <- g + geom_segment(aes(x = o$P_TRT_quantiles["50%", "PT_low"], y=0, 
                                  xend = o$P_TRT_quantiles["50%", "PT_low"], yend=1.05), linetype = 4)
        names(g$layers)[length(g$layers)] <- "low.PT.median"
        g <- g + geom_segment(aes(x = o$P_TRT_quantiles["50%", "lower.limit.TRT_low"], y=0, 
                                  xend = o$P_TRT_quantiles["50%", "lower.limit.TRT_low"], yend=1.15), linetype = 3)
        names(g$layers)[length(g$layers)] <- "low.TRT.limit.median"
        g <- g + geom_segment(aes(x = o$P_TRT_quantiles["50%", "higher.limit.TRT_low"], y=0, 
                                  xend = o$P_TRT_quantiles["50%", "higher.limit.TRT_low"], yend=1.15), linetype = 3)
        names(g$layers)[length(g$layers)] <- "high.TRT.limit.median"
        g <- g + geom_text(aes(x=o$P_TRT_quantiles["50%", "PT_low"], y=1.08), label = lab.PT)
        names(g$layers)[length(g$layers)] <- "low.PT.label"
        g <- g + geom_text(aes(x=o$P_TRT_quantiles["50%", "PT_low"], y=1.15), label = lab.TRT)
        names(g$layers)[length(g$layers)] <- "low.TRT.label"
        
        g <- g + geom_rect(aes(xmin = o$P_TRT_quantiles["50%", "lower.limit.TRT_high"], 
                               ymin = 0, 
                               xmax= o$P_TRT_quantiles["50%", "higher.limit.TRT_high"], 
                               ymax = 1), colour = NA, fill = col.TRT)
        names(g$layers)[length(g$layers)] <- "high.TRT.median"
        # limites de la limite basse de la TRT
        g <- g + geom_rect(aes(xmin = o$P_TRT_quantiles["2.5%", "lower.limit.TRT_high"], 
                               ymin = 0, 
                               xmax= o$P_TRT_quantiles["97.5%", "lower.limit.TRT_high"], 
                               ymax = 1), colour = NA, fill = col.TRT.CI)
        names(g$layers)[length(g$layers)] <- "high.TRT.lower.limit.credible.interval"
        # limites de la limite haute de la TRT
        g <- g + geom_rect(aes(xmin = o$P_TRT_quantiles["2.5%", "higher.limit.TRT_high"], 
                               ymin = 0, 
                               xmax= o$P_TRT_quantiles["97.5%", "higher.limit.TRT_high"], 
                               ymax = 1), colour = NA, fill = col.TRT.CI)
        names(g$layers)[length(g$layers)] <- "high.TRT.higher.limit.credible.interval"
        # limites de la PT
        g <- g + geom_rect(aes(xmin = o$P_TRT_quantiles["2.5%", "PT_high"], 
                               ymin = 0, 
                               xmax= o$P_TRT_quantiles["97.5%", "PT_high"], 
                               ymax = 1), colour = NA, fill = col.PT.CI)
        names(g$layers)[length(g$layers)] <- "high.PT.credible.interval"
        g <- g + geom_segment(aes(x = o$P_TRT_quantiles["50%", "PT_high"], y=0, 
                                  xend = o$P_TRT_quantiles["50%", "PT_high"], yend=1.05), linetype = 4)
        names(g$layers)[length(g$layers)] <- "high.PT.median"
        g <- g + geom_segment(aes(x = o$P_TRT_quantiles["50%", "lower.limit.TRT_high"], y=0, 
                                  xend = o$P_TRT_quantiles["50%", "lower.limit.TRT_high"], yend=1.15), linetype = 3)
        names(g$layers)[length(g$layers)] <- "high.TRT.lower.limit.median"
        g <- g + geom_segment(aes(x = o$P_TRT_quantiles["50%", "higher.limit.TRT_high"], y=0, 
                                  xend = o$P_TRT_quantiles["50%", "higher.limit.TRT_high"], yend=1.15), linetype = 3)
        names(g$layers)[length(g$layers)] <- "high.TRT.higher.limit.median"
        g <- g + geom_text(aes(x=o$P_TRT_quantiles["50%", "PT_high"], y=1.08), label = lab.PT)
        names(g$layers)[length(g$layers)] <- "PT.high.label"
        g <- g + geom_text(aes(x=o$P_TRT_quantiles["50%", "PT_high"], y=1.15), label = lab.TRT)
        names(g$layers)[length(g$layers)] <- "TRT.high.label"
        
      }
    }
  }
  if (show.observations) {
    if (show.observations.sd) {
    if (all(names(c(x$par, x$fixed.parameters)) != "n")) {
      if (males.freq) {  
        b <- getFromNamespace(".BinomialConfidence", ns="HelpersMG")(males,N)
        L1 <- modifyList(list(x=temperatures, y=males/N, bty="n", type="p", ylim=c(0,1), y.plus = b[,3], y.minus = b[,2]), L)
      } else {
        b <- getFromNamespace(".BinomialConfidence", ns="HelpersMG")(females,N)
        L1 <- modifyList(list(x=temperatures, y=females/N, bty="n", type="p", ylim=c(0,1), y.plus = b[,3], y.minus = b[,2]), L)
      }
    } else {
      if (males.freq) {  
        L1 <- modifyList(list(x=temperatures, y=males/N, bty="n", type="p", ylim=c(0,1)), L)
      } else {
        L1 <- modifyList(list(x=temperatures, y=females/N, bty="n", type="p", ylim=c(0,1)), L)
      }
    }
    L1 <- modifyList(L1, list(ylim=c(0,1), xlab="", ylab="", 
                              main="", axes=FALSE, xlim=c(x1, x2)))
    
    if (!use.ggplot) {
      par(xpd=FALSE)
      par(new=TRUE)
      
      a <- do.call(plot_errbar, L1) 
    } else {
      pch <- ifelse(is.null(L1$pch), 19, L1$pch)
      cex <- ifelse(is.null(L1$cex), 3, L1$cex)
      vx <- L1$x
      vy <- L1$y
      vymin <- L1$y.minus
      vymax <- L1$y.plus
      g <- g + geom_point(aes(x = vx, y=vy), shape = pch, size = cex)
      names(g$layers)[length(g$layers)] <- "Observations"
      g <- g + geom_errorbar(aes(x = vx, ymin = vymin, ymax = vymax))
      names(g$layers)[length(g$layers)] <- "Uncertainty.bars.observations"
    }
  }
  else {
    # J'affiche les points mais pas les observations
      if (males.freq) {  
        L1 <- modifyList(list(x=temperatures, y=males/N, bty="n", type="p", ylim=c(0,1)), L)
      } else {
        L1 <- modifyList(list(x=temperatures, y=females/N, bty="n", type="p", ylim=c(0,1)), L)
      }

    L1 <- modifyList(L1, list(ylim=c(0,1), xlab="", ylab="", 
                              main="", axes=FALSE, xlim=c(x1, x2)))
    
    if (!use.ggplot) {
      par(xpd=FALSE)
      par(new=TRUE)
      
      a <- do.call(plot_errbar, L1) 
    } else {
      # Pas sûr. Ca va peut être faire une erreur
      pch <- ifelse(is.null(L1$pch), 19, L1$pch)
      cex <- ifelse(is.null(L1$cex), 3, L1$cex)
      vx <- L1$x
      vy <- L1$y
      vymin <- L1$y.minus
      vymax <- L1$y.plus
      g <- g + geom_point(aes(x = vx, y=vy), shape = pch, size = cex)
      names(g$layers)[length(g$layers)] <- "Observations"
    }
  }
  }
  
  out_sr <- NULL
  
  if ((!is.null(o$sexratio_quantiles)) & (show.model)) {
    
    out_sr <- t(o$sexratio_quantiles)
    out_sr <- cbind(Temperatures=as.numeric(colnames(o$sexratio_quantiles)), 
                    out_sr)
    
    xi <- as.numeric(colnames(o$sexratio_quantiles))
    p <- o$sexratio_quantiles[2, ]
    
    
    if (males.freq) {   
      L1 <- modifyList(list(x=xi, y=p, bty="n"), L)
    } else {
      L1 <- modifyList(list(x=xi, y=1-p, bty="n"), L)
    }
    L1 <- modifyList(L1, list(ylim=c(0,1), axes=FALSE, xlab="", ylab="", type="l", main="", xlim=c(x1, x2)))
    
    
    if (!use.ggplot) {
      L2 <- L1[!(names(L1) %in% c("errbar.tick", "errbar.lwd"
                                  , "errbar.lty", "errbar.col"
                                  , "errbar.y.polygon"
                                  , "errbar.y.polygon.list"))]
      par(new=TRUE)
      par(xpd=FALSE)
      a <- do.call(plot, L2)
    } else {
      vx1 <- L1$x
      vy1 <- L1$y
      g <- g + geom_line(aes(x = vx1, y=vy1), linetype = 1, inherit.aes = FALSE)
      names(g$layers)[length(g$layers)] <- "model.median"
    }
    
    if (show.CI) {
      
      pm <- o$sexratio_quantiles[1, ]
      pp <- o$sexratio_quantiles[3, ]
      
      
      if (males.freq) {
        L1 <- modifyList(list(x=xi, y=pm, bty="n"), L)
      } else {
        L1 <- modifyList(list(x=xi, y=1-pm, bty="n"), L)
      }
      L1 <- modifyList(L1, list(ylim=c(0,1), axes=FALSE, xlab="", ylab="", type="l", main="", lty=2, xlim=c(x1, x2)))
      
      if (!use.ggplot) {
        L2 <- L1[!(names(L1) %in% c("errbar.tick", "errbar.lwd"
                                    , "errbar.lty", "errbar.col"
                                    , "errbar.y.polygon"
                                    , "errbar.y.polygon.list"))]
        par(new=TRUE)
        par(xpd=FALSE)
        a <- do.call(plot, L2) 
      } else {
        vx2 <- L1$x
        vy2 <- L1$y
        g <- g + geom_line(aes(x = vx2, y=vy2), linetype = L1$lty, inherit.aes = FALSE)
        names(g$layers)[length(g$layers)] <- "model.lower.credible.interval"
      }
      
      if (males.freq) {   
        L1 <- modifyList(list(x=xi, y=pp, bty="n"), L)
      } else {
        L1 <- modifyList(list(x=xi, y=1-pp, bty="n"), L)
      }
      L1 <- modifyList(L1, list(ylim=c(0,1), axes=FALSE, xlab="", ylab="", type="l", main="", lty=2, xlim=c(x1, x2)))
      
      if (!use.ggplot) {
        L2 <- L1[!(names(L1) %in% c("errbar.tick", "errbar.lwd"
                                    , "errbar.lty", "errbar.col"
                                    , "errbar.y.polygon"
                                    , "errbar.y.polygon.list"))]
        par(new=TRUE)
        par(xpd=FALSE)
        a <- do.call(plot, L2) 
      } else {
        vx3 <- L1$x
        vy3 <- L1$y
        g <- g + geom_line(aes(x = vx3, y=vy3), linetype = L1$lty, inherit.aes = FALSE)
        names(g$layers)[length(g$layers)] <- "model.upper.credible.interval"
      }
      
    }
  }
  
  if (use.ggplot) {
    plot(g)
    return(g)
  } else {
    
    return(invisible(out_sr))
  }
}





