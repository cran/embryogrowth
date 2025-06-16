#' HeterogeneityNests models heterogeneity of temperatures.
#' @title Model heterogeneity of temperatures.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Nothing.
#' @param nests An object of class Nests, Nests2, NestsResult, or mcmcComposite.
#' @param control.legend.total A list of options for legend.
#' @param control.legend.metabolicheating A list of options for legend.
#' @param show.full.incubation Show the plot with full incubation?
#' @param show.first.half.incubation Show the plot with first half incubation?
#' @param ... Parameters used for plot.
#' @description Generate a model of heterogeneity of temperatures.\cr
#' @examples 
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
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
#' nests <- FormatNests(data=nest, previous=NULL, col.Time="Time", LayingTime=Laying.Time_f)
#' HeterogeneityNests(nests, ylim=c(0, 4))
#' }
#' @export


HeterogeneityNests <-
  function(nests=stop("An object of class Nests2, NestsResult, or mcmcComposite."), 
           control.legend.total=list(), 
           control.legend.metabolicheating=list(), 
           show.full.incubation = TRUE, 
           show.first.half.incubation = TRUE, 
           ...) {
    
    p3p <- tryCatch(list(...), error=function(e) list())
    
    if (inherits(nests, "mcmcComposite")) {
      nests <- nests$parametersMCMC$control$temperatures
    }
    
    if (inherits(nests, "Nests")) {
      nests <- UpdateNests(nests)
    }
    
    if (inherits(nests, "NestsResult")) {
      nests <- nests$data
    }
    
    # Model
    modelH <- function(par, data, fixed.parameters=NULL) {
      par <- c(par, fixed.parameters)
      theoric <- par["min"]+(par["max"]-par["min"])*(1/(1+exp((1/par["S"])*(par["P"]-data[ ,"nb"]))))
      if (any(abs(par["asd"])*theoric+abs(par["bsd"]) < 0)) return(1E9)
      return(sum(-dnorm(x=data[ ,"H"], mean=theoric, sd=abs(par["asd"])*theoric+abs(par["bsd"]), log = TRUE)))
    }
    
    name_nests <- nests$Names
    
    LayingTime <- lapply(nests$Nests, FUN=function(x) unname(as.POSIXct(x$LayingTime)))
    LayingTime <- do.call("c", LayingTime)
    
    if (is.null(LayingTime)) {
      stop("Heterogeneity.Nests requires LayingTime being defined in FormatNests()")
    }
    
    BeginTime <- as.POSIXlt(min(LayingTime))
    BeginTime$mday <- 1
    BeginTime$mon <- 0
    BeginTime$hour <- 0
    BeginTime$min <- 0
    BeginTime$sec <- 0
    
    EndTime <- lapply(nests$Nests, FUN=function(x) unname(as.POSIXct(x$LayingTime) + max(x$data[, "Time"])*60))
    EndTime <- do.call("c", EndTime)
    
    EndTime <- as.POSIXlt(max(EndTime))
    EndTime$mday <- 1
    EndTime$mon <- 0
    EndTime$hour <- 0
    EndTime$min <- 0
    EndTime$sec <- 0
    
    
    dates <- lapply(name_nests, FUN = function(x) as.POSIXct(nests$Nests[[x]]$data[, "Time"] * 60 + LayingTime[x]))
    dates <- do.call("c", dates)
    
    dates <- sort(unique(dates))
    tempsmatrix <- matrix(data=NA, nrow=length(name_nests), ncol = length(dates))
    
    BeginTime <- as.POSIXct(BeginTime)
    EndTime <- as.POSIXct(EndTime)
    
    for (j in name_nests) {
      Eret <- nests$Nests[[j]]$data
      d <- nests$Nests[[j]]$LayingTime + Eret[, "Time"]*60
      tempsmatrix[which(j == name_nests), match(d, dates)] <- Eret[, "Temperatures C"]
    }
    
    # mean_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) mean(x, na.rm=TRUE))
    # sd_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) sd(x, na.rm=TRUE))
    qt_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) quantile(x, probs=c(0.025, 0.5, 0.975), na.rm=TRUE))
    qt_median_ma <- MovingWindow(qt_temp["50%", , drop=TRUE], window = 10, hole = "both")
    qt_low_ma <- MovingWindow(qt_temp["2.5%", , drop=TRUE], window = 10, hole = "both")
    qt_high_ma <- MovingWindow(qt_temp["97.5%", , drop=TRUE], window = 10, hole = "both")
    
    qt_median_ma <- qt_temp["50%", , drop=TRUE]
    qt_low_ma <- qt_temp["2.5%", , drop=TRUE]
    qt_high_ma <- qt_temp["97.5%", , drop=TRUE]
    
    
    nb_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) length(x[!is.na(x)]))
    
    dta <- data.frame(nb=nb_temp, H=qt_high_ma-qt_low_ma)
    dta <- dta[dta$nb > 1, ]
    
    mH <- optim(par=c('min' = 1.4111407511645899,
                      'max' = 3.9438420631957714,
                      'S' = 0.62127959812443434,
                      'P' = 5.6986305615475654,
                      'asd' = 0.13984929749918499,
                      'bsd' = 1.1902517432762878),
                fn=modelH, data=dta, hessian = TRUE,
                fixed.parameters=NULL,
                control = list(maxit=5000))
    par <- mH$par
    
    # MCMC estimation of parameters
    priors <- setPriors(
      par = par,
      se = SEfromHessian(mH$hessian),
      density = "dunif",
      rules = rules <- rbind(data.frame(Name="min", Min=0, Max=1),
                             data.frame(Name="max", Min=0, Max=10),
                             data.frame(Name="S", Min=0, Max=10),
                             data.frame(Name="P", Min=0, Max=20),
                             data.frame(Name="asd", Min=0, Max=2),
                             data.frame(Name="bsd", Min=0, Max=1)),
      silent = TRUE
    )
    
    priors$SDProp <- c('min' = 0.0041631916553423369,
                       'max' = 0.068167756161994636,
                       'S' = 0.009489156869831249,
                       'P' = 0.015037055049696076,
                       'asd' = 0.0078597369537145526,
                       'bsd' = 0.00023341615085828698)
    
    mcmc <- MHalgoGen(likelihood = modelH, parameters=priors, parameters_name = "par",
                      n.adapt = 0, thin = 5, n.iter = 50000,
                      data=dta, fixed.parameter=NULL,  trace = FALSE, adaptive = FALSE)
    
    se <- SEfromHessian(mH$hessian)
    
    x <- seq(from=2, to=max(nb_temp), by=0.1)
    
    outpar <- RandomFromHessianOrMCMC(
      mcmc = mcmc,
      chain = 1,
      regularThin = TRUE,
      MinMax = NULL,
      fitted.parameters = NULL,
      fixed.parameters = NULL,
      method = "MCMC",
      replicates = 10000,
      probs = c(0.025, 0.5, 0.975),
      fn = function(par) {y <- (par["min"]+(par["max"]-par["min"])*(1/(1+exp((1/par["S"])*(par["P"]-x)))));
      return(rnorm(length(y), y,  abs(par["asd"]) *y + abs(par["bsd"]) ))},
      silent = TRUE,
      ParTofn = "par"
    )
    
    outpar_se <- RandomFromHessianOrMCMC(
      mcmc = mcmc,
      chain = 1,
      regularThin = TRUE,
      MinMax = NULL,
      fitted.parameters = NULL,
      fixed.parameters = NULL,
      method = "MCMC",
      replicates = 10000,
      probs = c(0.025, 0.5, 0.975),
      fn = function(par) {y <- (par["min"]+(par["max"]-par["min"])*(1/(1+exp((1/par["S"])*(par["P"]-x)))));
      return(y)},
      silent = TRUE,
      ParTofn = "par"
    )
    
    mcmc_tot <- mcmc
    
    par(mar=c(4, 4, 1, 1))
    
    param.plot <- modifyList(modifyList(list(x=dta$nb, y=dta$H, xlim=c(2, max(dta$nb)), ylim=c(0, 8),
                                             pch=19, xlab="Number of timeseries", ylab = expression("Heterogeneity in "*~degree*C*""),
                                             bty="n", las=1, col=rgb(blue=0.1, green=0.1, red=1, alpha=0.01), xaxt="n"), 
                                        p3p), list(type = "n"))
    
    do.call("plot", param.plot)
    axis(1, at=2:max(nb_temp))
    
    if (show.full.incubation) {
      
      points(x=dta$nb, y=dta$H, col=rgb(blue=0.1, green=0.1, red=1, alpha=0.01), pch=19)
      
      # lines(x, outpar$quantiles["50%", ])
      polygon(x=c(x, rev(x)), y=c(outpar$quantiles["2.5%", ], rev(outpar$quantiles["97.5%", ])),
              border = NA, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.2))
      polygon(x=c(x, rev(x)), y=c(outpar_se$quantiles["2.5%", ], rev(outpar_se$quantiles["97.5%", ])),
              border = NA, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.6))
      
      segments(x0=1, x1=max(nb_temp), y0=summary(mcmc)$quantiles["max", "50%"], y1=summary(mcmc)$quantiles["max", "50%"], col="red", lty=4, lwd=2)
      segments(x0=1, x1=max(nb_temp), y0=summary(mcmc)$quantiles["max", "2.5%"], y1=summary(mcmc)$quantiles["max", "2.5%"], col="red", lty=3, lwd=2)
      segments(x0=1, x1=max(nb_temp), y0=summary(mcmc)$quantiles["max", "97.5%"], y1=summary(mcmc)$quantiles["max", "97.5%"], col="red", lty=3, lwd=2)
      cat("Taking into account all data in time series\n")
      cat("-------------------------------------------\n")
      cat(paste("The model converges to an heterogeneity of temperatures of",
                specify_decimal(unname(summary(mcmc)$quantiles["max", "50%"]), decimals=2),
                "\u00B0C (95% Credible Interval between",
                specify_decimal(unname(summary(mcmc)$quantiles["max", "2.5%"]), decimals=2),
                "and", specify_decimal(unname(summary(mcmc)$quantiles["max", "97.5%"]), decimals=2),
                ") at the scale of the beach.\n"))
      cat(paste("The maximum heterogeneity is",
                specify_decimal(unname(summary(mcmc)$quantiles["bsd", "50%"])+unname(summary(mcmc)$quantiles["asd", "50%"])*unname(summary(mcmc)$quantiles["max", "50%"]) * 1.96 + unname(summary(mcmc)$quantiles["max", "50%"]), decimals=2),
                "\u00B0C.\n"))
      
    }
    
    if (show.first.half.incubation) {
      
      # With only first half of timeseries
      tempsmatrix[] <- NA
      for (j in name_nests) {
        Eret <- nests$Nests[[j]]$data
        maxt <- max(Eret[, "Time"])/2
        Eret <- Eret[Eret[, "Time"] < maxt, ]
        d <- nests$Nests[[j]]$LayingTime + Eret[, "Time"]*60
        
        tempsmatrix[which(j == name_nests), match(d, dates)] <- Eret[, "Temperatures C"]
      }
      
      # mean_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) mean(x, na.rm=TRUE))
      # sd_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) sd(x, na.rm=TRUE))
      qt_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) quantile(x, probs=c(0.025, 0.5, 0.975), na.rm=TRUE))
      qt_median_ma <- MovingWindow(qt_temp["50%", , drop=TRUE], window = 10, hole = "both")
      qt_low_ma <- MovingWindow(qt_temp["2.5%", , drop=TRUE], window = 10, hole = "both")
      qt_high_ma <- MovingWindow(qt_temp["97.5%", , drop=TRUE], window = 10, hole = "both")
      
      qt_median_ma <- qt_temp["50%", , drop=TRUE]
      qt_low_ma <- qt_temp["2.5%", , drop=TRUE]
      qt_high_ma <- qt_temp["97.5%", , drop=TRUE]
      
      
      nb_temp <- apply(tempsmatrix, MARGIN = 2, FUN = function(x) length(x[!is.na(x)]))
      
      dta <- data.frame(nb=nb_temp, H=qt_high_ma-qt_low_ma)
      dta <- dta[dta$nb > 1, ]
      
      pari <- as.parameters(mcmc, index="median")
      
      
      mH <- optim(par=pari[c("min", "max", "asd", "bsd")],
                  fn=modelH, data=dta, hessian = TRUE,
                  fixed.parameters=pari[c("S", "P")],
                  control = list(maxit=5000))
      par <- mH$par
      
      # MCMC estimation of parameters
      priors <- setPriors(
        par = par,
        se = SEfromHessian(mH$hessian),
        density = "dunif",
        rules = rules <- rbind(data.frame(Name="min", Min=0, Max=1),
                               data.frame(Name="max", Min=0, Max=10),
                               data.frame(Name="S", Min=0, Max=10),
                               data.frame(Name="P", Min=0, Max=20),
                               data.frame(Name="asd", Min=0, Max=2),
                               data.frame(Name="bsd", Min=0, Max=1)),
        silent = TRUE
      )
      
      
      
      priors$SDProp <- c('min' = 0.0041631916553423369,
                         'max' = 0.068167756161994636,
                         'asd' = 0.0078597369537145526,
                         'bsd' = 0.00023341615085828698)
      mcmc <- MHalgoGen(likelihood = modelH, parameters=priors, parameters_name = "par",
                        n.adapt = 0, thin = 5, n.iter = 50000,
                        fixed.parameters=pari[c("S", "P")],
                        data=dta, trace = FALSE, adaptive = FALSE)
      
      se <- SEfromHessian(mH$hessian)
      
      x <- seq(from=2, to=max(nb_temp), by=0.1)
      
      outpar <- RandomFromHessianOrMCMC(
        mcmc = mcmc,
        chain = 1,
        regularThin = TRUE,
        MinMax = NULL,
        fitted.parameters = NULL,
        fixed.parameters = pari[c("S", "P")],
        method = "MCMC",
        replicates = 10000,
        probs = c(0.025, 0.5, 0.975),
        fn = function(par) {y <- (par["min"]+(par["max"]-par["min"])*(1/(1+exp((1/par["S"])*(par["P"]-x)))));
        return(rnorm(length(y), y, abs(par["asd"]) * y + abs(par["bsd"])))},
        silent = TRUE,
        ParTofn = "par"
      )
      
      outpar_se <- RandomFromHessianOrMCMC(
        mcmc = mcmc,
        chain = 1,
        regularThin = TRUE,
        MinMax = NULL,
        fitted.parameters = NULL,
        fixed.parameters = pari[c("S", "P")],
        method = "MCMC",
        replicates = 10000,
        probs = c(0.025, 0.5, 0.975),
        fn = function(par) {y <- (par["min"]+(par["max"]-par["min"])*(1/(1+exp((1/par["S"])*(par["P"]-x)))));
        return(y)},
        silent = TRUE,
        ParTofn = "par"
      )
      
      
      points(x=dta$nb, y=dta$H, pch=19,
             col=rgb(blue=1, green=0.1, red=0.1, alpha=0.01))
      
      polygon(x=c(x, rev(x)), y=c(outpar$quantiles["2.5%", ], rev(outpar$quantiles["97.5%", ])),
              border = NA, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2))
      polygon(x=c(x, rev(x)), y=c(outpar_se$quantiles["2.5%", ], rev(outpar_se$quantiles["97.5%", ])),
              border = NA, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.6))
      
      pari <- as.parameters(mcmc, index="median")
      
      segments(x0=2, x1=max(nb_temp), y0=unname(summary(mcmc)$quantiles["max", "50%"]), y1=unname(summary(mcmc)$quantiles["max", "50%"]), col="blue", lty=4, lwd=2)
      
      segments(x0=2, x1=max(nb_temp), y0=summary(mcmc)$quantiles["max", "2.5%"], y1=summary(mcmc)$quantiles["max", "2.5%"], col="blue", lty=3, lwd=2)
      segments(x0=2, x1=max(nb_temp), y0=summary(mcmc)$quantiles["max", "97.5%"], y1=summary(mcmc)$quantiles["max", "97.5%"], col="blue", lty=3, lwd=2)
      
      cat("Taking into account only first half time series\n")
      cat("-----------------------------------------------\n")
      cat(paste("The model converges to an heterogeneity of temperatures of",
                specify_decimal(unname(summary(mcmc)$quantiles["max", "50%"]), decimals=2),
                "\u00B0C (95% Credible Interval between",
                specify_decimal(unname(summary(mcmc)$quantiles["max", "2.5%"]), decimals=2),
                "and", specify_decimal(unname(summary(mcmc)$quantiles["max", "97.5%"]), decimals=2),
                ") at the scale of the beach.\n"))
      cat(paste("The maximum heterogeneity is",
                specify_decimal(unname(summary(mcmc)$quantiles["bsd", "50%"])+unname(summary(mcmc)$quantiles["asd", "50%"])*unname(summary(mcmc)$quantiles["max", "50%"]) * 1.96 + unname(summary(mcmc)$quantiles["max", "50%"]), decimals=2),
                "\u00B0C.\n"))
    }
    #
    control.legend.total <- modifyList(list(x="bottomright", legend = c("Observations", "95% SD", "95% Bayesian CI model", "Heterogeneity of the beach"),
                                            col=c("lightgrey", "grey", "black", "black"), 
                                            pch=c(19, 15, NA, NA), lty=c(NA, NA, 1, 3), lwd=c(NA, NA, 3, 1)), 
                                       control.legend.total)
    do.call("legend", control.legend.total)
    
    control.legend.metabolicheating <- modifyList(list(x="topleft", legend = c("With Metabolic heating", "Without metabolic heating"), 
                                                       lty=c(1, 1), col=c("red", "blue")), 
                                                  control.legend.metabolicheating)
    
    do.call("legend", control.legend.metabolicheating)
    
    return(invisible(c(lower.max.half=unname(summary(mcmc)$quantiles["max", "2.5%"]), 
                       median.max.half=unname(summary(mcmc)$quantiles["max", "50%"]), 
                       upper.max.half=unname(summary(mcmc)$quantiles["max", "97.5%"]), 
                       lower.max=unname(summary(mcmc_tot)$quantiles["max", "2.5%"]), 
                       median.max=unname(summary(mcmc_tot)$quantiles["max", "50%"]), 
                       upper.max=unname(summary(mcmc_tot)$quantiles["max", "97.5%"])
    )
    )
    )
  }
