#' plot.Nests2 shows the plot of temperatures for set of nests
#' @title Show the plot of temperatures with set of nests
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @return The position of labels if xaxt="n" is used.
#' @param x Data formated using formatdata.
#' @param series Series to be used, logical (TRUE ou FALSE), numbers or names. If "all", all series are used.
#' @param time Can be relative or absolute.
#' @param show.heterogeneity If TRUE, show the 95% heterogeneity in grey.
#' @param show.legend.heterogeneity Show the heterogeneity legend.
#' @param col.heterogeneity Color of heterogeneity.
#' @param probs.heterogeneity Quantiles of heterogeneity.
#' @param show.legend.names Show a legend for names.
#' @param control.legend.heterogeneity The list of parameters for legend of heterogeneity.
#' @param control.legend.names The list of parameters for legend of names.
#' @param col A recycled vector of colors.
#' @param xlim The xlim parameter of plot or can be "auto", "auto-month", or "auto-year"
#' @param xlab.axis Labels of x-axis
#' @param cex.axis The size of x-axis labels
#' @param ... Parameters used by plot function
#' @description Show the plot of temperatures with set of nests\cr
#' If time is "absolute", LayingTime must be indicated in FormatNests()\cr
#' xlim being auto-month used the closest begin and end of month and 
#' auto-year uses the clsest begin and end of year.
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' plot(x=formated, series="all", col=rainbow(21))
#' plot(x=formated, series=1, main="DY.1")
#' Laying.Time <- matrix(c("DY.1", "15/05/2010", 
#'                  "DY.17", "24/05/2010", 
#'                  "DY.16", "24/05/2010", 
#'                  "DY.18", "25/05/2010", 
#'                  "DY.20", "25/05/2010", 
#'                  "DY.21", "26/05/2010", 
#'                  "DY.22", "26/05/2010", 
#'                  "DY.23", "26/05/2010", 
#'                  "DY.24", "27/05/2010", 
#'                  "DY.25", "27/05/2010", 
#'                  "DY.28", "28/05/2010", 
#'                  "DY.26", "28/05/2010", 
#'                  "DY.27", "28/05/2010", 
#'                  "DY.146", "20/06/2010", 
#'                  "DY.147", "20/06/2010", 
#'                  "DY.172", "24/06/2010", 
#'                  "DY.175", "24/06/2010", 
#'                  "DY.170", "24/06/2010", 
#'                  "DY.260", "06/07/2010", 
#'                  "DY.282", "12/07/2010", 
#'                  "DY.310", "18/07/2010", 
#'                  "DY.309", "18/07/2010", 
#'                  "DY.328", "25/07/2010", 
#'                  "DY.331", "26/07/2010"), byrow=TRUE, ncol=2)
#' tz <- OlsonNames()[grepl("Asia/Istanbul", OlsonNames())]
#' Laying.Time_f <- as.POSIXlt.character(Laying.Time[, 2], format = "%d/%m/%Y", tz=tz)
#' names(Laying.Time_f) <- Laying.Time[, 1]
#' formated <- FormatNests(data=nest, previous=NULL, col.Time="Time", LayingTime=Laying.Time_f)
#' plot(x=formated, time="absolute", ylim=c(20, 35), 
#'      col= rainbow(21, alpha = 1), control.legend.heterogeneity=list(cex.0.5))
#' }
#' @method plot Nests2
#' @export

plot.Nests2 <- function(x                                                                     , 
                        series="all"                                                          , 
                        time="relative"                                                       , 
                        show.heterogeneity=TRUE                                               , 
                        probs.heterogeneity=c(0.025, 0.5, 0.975)                              , 
                        col.heterogeneity=rgb(red = 0.2, green = 0.2, blue = 0.2, alpha = 0.6), 
                        show.legend.names=TRUE                                                , 
                        show.legend.heterogeneity=TRUE                                        , 
                        control.legend.heterogeneity=list()                                   , 
                        control.legend.names=list()                                           , 
                        col="black"                                                           ,
                        cex.axis=1                                                            ,
                        xlim="auto"                                                           ,
                        xlab.axis=NULL                                                        ,
                        ...                                                                   ) {
  
  p3p <- tryCatch(list(...), error=function(e) list())
  
  atT <- NULL
  
  if (is.character(xlim)) {
    xlim <- tolower(xlim)
    xlim <- match.arg(xlim, c("auto", "auto-month", "auto-year"))
    if (xlim == "auto") xlim <- "auto-month"
  }
  
  if ((show.heterogeneity) & (length(probs.heterogeneity) != 3) & (time == "absolute")) 
    stop("probs.heterogeneity must be of length 3.")
  
  time <- match.arg(tolower(time), choices = c("relative", "absolute"))
  
  names_nests <- x$Names
  
  if (series[[1]] == "all") {
    series <- names_nests
  } else {
    if (!is.character(series)) series <- names_nests[series]
  }
  
  
  if (length(series) == 0) {
    warning("No series to plot")
    return()
  } 
  
  col <- rep(col, length(series))[1:length(series)]
  
  if (length(p3p[[c("lty")]]) != 0) {
    lty <- rep(p3p[[c("lty")]], length(series))[1:length(series)]
    p3p <- p3p[names(p3p) != "lty"]
  } else {
    lty <- rep(1, length(series))
  }
  if (length(p3p[[c("lwd")]]) != 0) {
    lwd <- rep(p3p[[c("lwd")]], length(series))[1:length(series)]
    p3p <- p3p[names(p3p) != "lwd"]
  } else {
    lwd <- rep(1, length(series))
  }
  
  nids <- x$Nests[series]
  
  if (time == "relative") {
    
    if (is.character(xlim)) {
      xlim <- c(0,70)
    }
    
    L <- list(xlab="Days of incubation", ylab=expression("Temperature in " * degree * "C"), 
              main="", bty="n", las=1, type="n", xlim=xlim, ylim=c(20, 35))
    
    L2 <- c(list(x=0, y=0), p3p)
    
    L3 <- modifyList(x=L, val = L2)
    
    a <- do.call(plot, L3) 
    
    
    index <- 1
    for (j in series) {
      nidsx <- nids[[j]]$data
      x <- nidsx[, "Time"]/60/24
      y <- nidsx[, "Temperatures C"]
      L <- modifyList(list(x=x, y=y), list(col=col[index], lty=lty[index], lwd=lwd[index]))
      index <- index + 1
      a <- do.call(lines, L) 
    }
  } else {
    # Je suis en absolute
    
    if (is.character(xlim)) {
      
      BeginTime <- lapply(nids, FUN=function(x) unname(as.POSIXct(x$LayingTime)))
      BeginTime <- do.call("c", BeginTime)
      
      if (is.null(BeginTime)) {
        stop("Absolute time requires LayingTime being defined in FormatNests()")
      }
      
      # Attention, si en secondes, ce n'est pas bon. Corrigé 5/5/2025
      EndTime <- lapply(series, FUN = function(x) {
        if ((nids[[x]][["UnitTime"]]) == "minutes") r <- as.POSIXct(nids[[x]][["data"]][, "Time"] * 60 + nids[[x]][["LayingTime"]])
        if ((nids[[x]][["UnitTime"]]) == "seconds") r <- as.POSIXct(nids[[x]][["data"]][, "Time"] + nids[[x]][["LayingTime"]])
        if ((nids[[x]][["UnitTime"]]) == "hours") r <- as.POSIXct(nids[[x]][["data"]][, "Time"] * 60*60 + nids[[x]][["LayingTime"]])
        return(unname(max(r)))
      }
      )
      
      EndTime <- do.call("c", EndTime)
      
      if (xlim == "auto-year") {
        
        BeginTime <- as.POSIXlt(min(BeginTime))
        BeginTime$mday <- 1
        BeginTime$mon <- 0
        BeginTime$hour <- 0
        BeginTime$min <- 0
        BeginTime$sec <- 0
        
        EndTime <- as.POSIXlt(max(EndTime))
        EndTime$mday <- 1
        EndTime$mon <- 0
        EndTime$hour <- 0
        EndTime$min <- 0
        EndTime$sec <- 0
        
        EndTime$year <- EndTime$year + 1
        
        BeginTime <- as.POSIXct(BeginTime)
        EndTime <- as.POSIXct(EndTime)
        
      } else {
        BeginTime <- as.POSIXlt(min(BeginTime))
        BeginTime$mday <- 1
        BeginTime$hour <- 0
        BeginTime$min <- 0
        BeginTime$sec <- 0
        BeginTime <- as.POSIXct(BeginTime)
        
        
        EndTime <- as.POSIXlt(max(EndTime))
        EndTime$mday <- 1
        if (EndTime$mon == 1) {
          EndTime$mon <- 0
          EndTime$yea <- EndTime$yea + 1
        } else {
          EndTime$mon <- EndTime$mon + 1
        }
        EndTime$hour <- 0
        EndTime$min <- 0
        EndTime$sec <- 0
        EndTime <- as.POSIXct(EndTime)
      }
      
    } else {
      BeginTime <- as.POSIXct(xlim[1])
      EndTime <- as.POSIXct(xlim[2])
    }
    
    dates <- lapply(series, FUN = function(x) {
      if ((nids[[x]][["UnitTime"]]) == "minutes") r <- as.POSIXct(nids[[x]][["data"]][, "Time"] * 60 + nids[[x]][["LayingTime"]])
      if ((nids[[x]][["UnitTime"]]) == "seconds") r <- as.POSIXct(nids[[x]][["data"]][, "Time"] + nids[[x]][["LayingTime"]])
      if ((nids[[x]][["UnitTime"]]) == "hours") r <- as.POSIXct(nids[[x]][["data"]][, "Time"] * 60*60 + nids[[x]][["LayingTime"]])
      return(r)
    }
    )
    
    dates <- do.call("c", dates)
    
    dates <- sort(unique(dates))
    tempsmatrix <- matrix(data=NA, nrow=length(series), ncol = length(dates))
    
    
    L <- list(xlab="Months", ylab=expression("Temperature in " * degree * "C"), 
              main="", bty="n", las=1, type="n", xlim=c(BeginTime, EndTime), 
              ylim=c(20, 40), xaxt="n")
    
    L2 <- c(list(x=0, y=0), p3p)
    
    L3 <- modifyList(x=L, val = L2)
    
    a <- do.call(plot, L3) 
    
    BeginTime <- L3$xlim[1]
    EndTime <- L3$xlim[2]
    
    atT <- seq(from=BeginTime, to=EndTime, by="1 month")
    atC <- rep("", length(atT))
    atY <- unique(c(1, which((as.POSIXlt(atT)$mon == 0))))
    atC[atY] <- paste0(as.character(1+as.POSIXlt(atT[atY])$mon), "-", as.character(1900+as.POSIXlt(atT[atY])$year))
    atC[which(atC=="")] <- as.character(1+as.POSIXlt(atT[which(atC=="")])$mon)
    
    # Si p3p est null
    p3p_xaxt <- ifelse(identical(p3p, list()), TRUE, FALSE)
    p3p_xaxt2 <- ifelse(p3p_xaxt, TRUE, ifelse(is.null(p3p$xaxt), TRUE, ifelse(p3p$xaxt == "n", FALSE, TRUE)))
    
    if (p3p_xaxt2) {
      axis(1, at=seq(from=BeginTime, to=EndTime, by="1 month"), labels = NA)
      axis(1, at=atT[atY], labels = NA, lwd.ticks = 3)
      if (!is.null(xlab.axis)) {
        axis(1, at=atT, labels = xlab.axis, lwd.ticks = NA, gap.axis = -1, cex.axis=cex.axis)
      } else {
        axis(1, at=atT, labels = atC, lwd.ticks = NA, gap.axis = -1, cex.axis=cex.axis)
      }
      
    }
    
    index <- 1
    
    for (c in series) {
      Eret <- nids[[c]]$data
      d <- nids[[c]]$LayingTime + Eret[, "Time"]*60
      lines(x=d, y=Eret[, "Temperatures C"], col=col[index])
      index <- index + 1
      tempsmatrix[which(c == series), match(d, dates)] <- Eret[, "Temperatures C"]
    }
    
    if (show.heterogeneity) {
      
      mean_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) mean(x, na.rm=TRUE))
      sd_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) sd(x, na.rm=TRUE))
      qt_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) quantile(x, probs=probs.heterogeneity, na.rm=TRUE))
      qt_median_ma <- MovingWindow(qt_temp[2, , drop=TRUE], window = 10, hole = "both")
      qt_low_ma <- MovingWindow(qt_temp[1, , drop=TRUE], window = 10, hole = "both")
      qt_high_ma <- MovingWindow(qt_temp[3, , drop=TRUE], window = 10, hole = "both")
      
      qt_median_ma <- qt_temp[2, , drop=TRUE]
      qt_low_ma <- qt_temp[1, , drop=TRUE]
      qt_high_ma <- qt_temp[3, , drop=TRUE]
      
      
      nb_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) length(x[!is.na(x)]))
      
      plot_errbar(x=dates, 
                  y=qt_median_ma, 
                  y.minus=qt_low_ma, 
                  y.plus=qt_high_ma, 
                  errbar.y.polygon =TRUE, type="n", 
                  errbar.y.polygon.list = list(col=col.heterogeneity, 
                                               border=NA), 
                  add = TRUE)
      if (show.legend.heterogeneity) {
        control.legend.heterogeneity <- modifyList(list(x="topleft", legend = "Quantiles", col=col.heterogeneity, 
                                                        pch=15, cex = 0.8), control.legend.heterogeneity)
        do.call("legend", control.legend.heterogeneity)
      }
    }
  }
  if (show.legend.names) {
    control.legend.names <- modifyList(list(x="bottomright", legend = series, 
                                            lty=lty, lwd=lwd, col=col), control.legend.names)
    do.call("legend", control.legend.names)
  }
  
  return(invisible(atT))
}


