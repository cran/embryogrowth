#' plotR shows the fitted growth rate dependent on temperature and the density of the mcmc
#' @title Show the fitted growth rate dependent on temperature and its density
#' @author Marc Girondot
#' @return The value of scaleY to be used with other plotR function
#' @param result A result object or a list of result objects
#' @param resultmcmc A result object from GRTN_MHmcmc() function
#' @param parameters A set of parameters - Has the priority over result
#' @param fixed.parameters A set of fixed parameters
#' @param temperatures A set of temperatures - Has the priority over result
#' @param hessian An hessian matrix
#' @param replicate.CI Number of replicates to estimate confidence interval with Hessian if delta method failed
#' @param curves What curves to show: "MCMC quantiles" or "MCMC mean-SD" based on mcmc or "ML" or "ML quantiles" for maximum-likelihood
#' @param colramp Ramp function accepting an integer as an argument and returning n colors.
#' @param bandwidth numeric vector (length 1 or 2) of smoothing bandwidth(s). If missing, a more or less useful default is used. bandwidth is subsequently passed to function bkde2D.
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param cex.lab cex value for axis
#' @param cex.axis cex value for axis
#' @param bty Box around the pot
#' @param las Orientation for labels in y axis
#' @param main Title of the graph
#' @param set.par 1 or 2 to designate with set of parameters to show
#' @param pch Character for outlayers
#' @param col The color of the lines
#' @param col.polygon The color of the polygon
#' @param polygon If TRUE, confidence interval is shown as a polygon with color
#' @param lty The type of lines
#' @param ltyCI The type of lines
#' @param lwd The type of lines
#' @param lwdCI The type of lines
#' @param ylim Range of values for y-axis
#' @param xlim Range of values for x-axis
#' @param by.temperature Step to built the temperatures
#' @param scaleY Scaling factor for y axis or "auto"
#' @param show.density TRUE or FALSE for use with Hessian or MCMC
#' @param probs Confidence or credibility interval to show
#' @param new Should the graphics be a new one (TRUE) or superimposed to a previous one (FALSE) 
#' @param show.hist TRUE or FALSE
#' @param ylimH Scale of histogram using ylimH=c(min, max)
#' @param atH Position of ticks for scale of histogram
#' @param ylabH Label for histogram scale
#' @param breaks See ?hist
#' @param log.hist SHould the y scale for hist is log ?
#' @param mar The value of par("mar"). If null, it will use default depending on show.dist. If NA, does not change par("mar").
#' @description Show the fitted growth rate dependent on temperature and its density.\cr
#' The curve "ML quantiles" is based on delta method.\cr
#' The curve "ML" just shows the fitted model.\cr
#' The curve "MCMC quantiles" uses the mcmc replicates to build the quantiles.\cr
#' The curve "MCMC mean-SD" uses the mcmc replicates to build a symetric credibility interval.\cr
#' The parameter curves is case insensitive. If only parameters is given, curves must be ML.
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' plotR(result = resultNest_4p_SSM4p, 
#'              resultmcmc=resultNest_mcmc_4p_SSM4p, 
#'              curves = "MCMC quantiles")
#' #################
#' plotR(resultmcmc=resultNest_mcmc_4p_SSM4p, 
#'              curves = "MCMC quantiles", show.density=TRUE)
#' #################
#' plotR(resultmcmc=resultNest_mcmc_4p_SSM4p, 
#'              curves = "MCMC quantiles", polygon=TRUE)
#' #################
#' plotR(resultmcmc=resultNest_mcmc_6p_SSM6p, ylim=c(0,4), 
#'       curves = "MCMC quantiles", polygon=TRUE, col.polygon = rgb(0, 1, 0, 1))
#' plotR(resultmcmc=resultNest_mcmc_4p_SSM4p,
#'        curves = "MCMC quantiles", polygon=TRUE, col.polygon = rgb(1, 0, 0, 0.5), new=FALSE)
#' legend("topleft", legend=c("SSM 4 parameters", "SSM 6 parameters"), 
#'         pch=c(15, 15), col=c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 1)))
#' #################
#' sy <- plotR(resultmcmc=resultNest_mcmc_4p_SSM4p, 
#'              curves = "MCMC quantiles", show.density=FALSE)
#' plotR(resultmcmc=resultNest_mcmc_6p_SSM6p, col="red",
#'              curves = "MCMC quantiles", show.density=FALSE, 
#'              new=FALSE, scaleY=sy)
#' #################
#' sy <- plotR(result=resultNest_6p_SSM6p, curves="ML", 
#'              show.hist = TRUE, new = TRUE)
#' plotR(result=resultNest_4p_SSM4p, curves="ML", scaleY=sy, 
#'              show.hist = FALSE, new = FALSE, col="red")
#' #################
#' plotR(result=resultNest_6p_SSM6p, curves="ML", 
#'              show.hist = TRUE, ylimH=c(0,1), atH=c(0, 0.1, 0.2))
#' ################
#' plotR(result = resultNest_4p_SSM4p, 
#'              resultmcmc=resultNest_mcmc_4p_SSM4p, 
#'              show.density = TRUE, 
#'              curves = "MCMC quantiles")
#' #################
#' plotR(result=resultNest_4p_SSM4p, 
#'              ylim=c(0, 4), curves="ML quantiles", scaleY=1E5)
#' #################             
#' plotR(result=resultNest_4p_SSM4p, show.hist = TRUE,
#'              ylim=c(0, 4), curves="ML quantiles", scaleY=1E5)
#' #################
#' plotR(resultmcmc=resultNest_mcmc_4p_SSM4p, 
#'              ylim=c(0, 4), curves = "MCMC quantiles", show.density=TRUE, scaleY=1E5)
#' }
#' @export


plotR <-
  function(result = NULL, 
           resultmcmc = NULL, 
           parameters = NULL, fixed.parameters = NULL, 
           temperatures  = NULL,
           curves = "ML quantiles", set.par=1, ylim=c(0, 5), xlim=c(20,35), 
           hessian = NULL, 
           replicate.CI = 1000, 
           cex.lab = par("cex"), cex.axis = par("cex"),
           scaleY="auto", lty=1, ltyCI=3, lwd=1, lwdCI=1, col = "black", 
           col.polygon="grey", polygon=FALSE, probs=0.95,
           colramp=colorRampPalette(c("white", rgb(red = 0.5, green = 0.5, blue = 0.5))), 
           bandwidth = c(0.1, 0.01), pch = "", main="",
           xlab = expression("Temperature in "*degree*"C"), 
           ylab = NULL, bty = "n", las = 1, 
           by.temperature=0.1, show.density=FALSE, 
           new=TRUE, show.hist=FALSE, ylimH = NULL, atH = NULL, 
           ylabH="Temperature density", breaks = "Sturges", 
           log.hist=FALSE,
           mar=NULL) {
    
    # result = NULL; resultmcmc=NULL; parameters = NULL; fixed.parameters = NULL; hessian = NULL; replicate.CI=1000; probs=0.95; temperatures  = NULL;curves = "ML quantiles"; set.par=1; ylim=c(0, 5); xlim=c(20,35); cex.lab = 1; cex.axis = 1;scaleY="auto"; lty=1; ltyCI=3; lwd=1; lwdCI=1; colramp=colorRampPalette(c("white", rgb(red = 0.5, green = 0.5, blue = 0.5))); bandwidth = c(0.3, 0.05); pch = ""; main=""; col = "black"; col.polygon="grey"; polygon=FALSE; xlab = expression("Temperature in"*degree*"C"); ylab = NULL; bty = "n"; las = 1; by.temperature=0.1; show.density=FALSE;new=TRUE;show.hist=FALSE; ylimH = NULL; atH = NULL;ylabH="Temperature density";breaks = "Sturges";log.hist=FALSE;mar=NULL
    # resultmcmc = resultNest_mcmc_4p_SSM4p
    # result = resultNest_4p_SSM4p
    
    curves <- tolower(curves)
    
    SSM <- getFromNamespace(".SSM", ns="embryogrowth")
    
    if (class(result)=="mcmcComposite") {
      resultmcmc <- result
      result <- NULL
    }
    
    if (class(result)=="numeric") {
      parameters <- result
      result <- NULL
    }
    
    if (is.null(hessian) & !is.null(result)) hessian <- result$hessian
    if (is.null(fixed.parameters) & !is.null(result)) fixed.parameters <- result$fixed.parameters
    if (is.null(parameters) & !is.null(result)) parameters <- result$par
    
#    parameters <- c(parameters, fixed.parameters)
#    if (is.null(parameters) & !is.null(result)) parameters <- c(result$par, result$fixed.parameters)
    if (is.null(parameters) & !is.null(resultmcmc)) parameters <- suppressMessages(as.parameters(resultmcmc))
    
    if (is.null(temperatures) & !is.null(result)) temperatures <- result$data
    
    if (show.density & is.null(resultmcmc) & (is.null(hessian))) {
      warning("show.density option needs or a mcmc or an Hessian object")
      show.density <- FALSE
    }
    
    if (show.hist & is.null(temperatures)) {
      warning("show.hist option needs a result object or temperatures")
      show.hist <- FALSE
    }
    
      if (is.null(mar)) {
        if (!is.null(result) & show.hist & new) {
          par(mar=c(4, 4, 1, 4)+0.4) 
        } else {
          par(mar=c(4, 4, 4, 1)+0.4)
        }
      } else {
        if (all(!is.na(mar))) par(mar=mar)
      }
    
    temp <- seq(from=xlim[1], to=xlim[2], by=by.temperature)
    
    rxxw_ML <- NULL
    rxxw <- NULL
    xxw <- NULL
    xxw_ML <- NULL
    
    if (!is.null(resultmcmc) & (show.density | curves=="mcmc quantiles" | curves=="mcmc mean-sd")) {
    dataR <- resultmcmc$resultMCMC[[1]]
    xxw <- matrix(ncol = 2, 
                  nrow = length(temp)*nrow(dataR))
    xxw2 <- matrix(data = as.numeric(), ncol = length(temp), 
                   nrow = nrow(dataR))
    
    xxw[, 1] <- rep(temp, nrow(dataR))
    for (i in 1:nrow(dataR)) {
      r <- SSM(T=temp, 
               parms=c(dataR[i,], fixed.parameters))[[set.par]]
      jhn <- (1+(i-1)*length(temp)):((i-1)*length(temp)+length(temp))
      # xxw[jhn, 1] <- temp
      xxw[jhn, 2] <- r
      xxw2[i, ] <- r
    }
    rxxw <- apply(xxw2, MARGIN = 2, 
                  FUN = function(x) c(Mean=mean(x), sd=sd(x), quantile(x, probs =c((1-probs)/2, 0.5, 1-(1-probs)/2), na.rm=TRUE)))
    rxxw <- rbind(temperatures=temp, rxxw)
    }
    
    if (curves=="ml quantiles" & is.null(hessian)) {
      curves <- "ml"
      warning("ml quantiles need an hessian matrix")
    }
    
    # j'ai un hessian qui est fourni
    if (curves=="ml quantiles" & !is.null(hessian)) {
      

      Scoefs <- parameters
#      hessian <- result$hessian
      # rownames(hessian) <- colnames(hessian) <- names(Scoefs)
      

      
      rxxw_ML <- data.frame(temperatures = numeric(), 
                            Mean=numeric(), sd=numeric(), 
                            "2.5%" = numeric(), 
                            "50%" = numeric(), 
                            "97.5%" = numeric()) 
      
      txt <- paste0("getFromNamespace('.wrapperSSM', ns='embryogrowth')(", 
      paste0("b[", 1:length(Scoefs), "], ", collapse=""), 
      ifelse(length(fixed.parameters)>0, 
      paste0("x[1]", ", x[", 2:(1+length(fixed.parameters)), "] ", collapse=""), 
      paste0("x[1]")), 
      ")", collapse="")
    warn <- FALSE
      for (temperature in temp) {
        options(warn=2) 
        # if it returns an error, it means that CI is too small to be estimated
        ic5 <- try(getFromNamespace(".nlConfint", ns="HelpersMG")(texts=txt, 
                                              level = probs, coeff = Scoefs,
                                              Vcov = solve(hessian), df2 = TRUE, 
                                              x = c(T=temperature, fixed.parameters), 
                                              silent=TRUE), silent=TRUE)
        options(warn=0)
        vy <- getFromNamespace('.wrapperSSM', ns='embryogrowth')(c(Scoefs, fixed.parameters, T=temperature))
        if (class(ic5) == "try-error") {
          # L'inervalle de confiance est nul ou autre probleme dans le genre
          # warning(ic5)
          warn <- TRUE
          # rxxw_ML <- rbind(rxxw_ML, 
          #                data.frame(temperatures=temperature, 
          #                           Mean = vy, sd = NA, 
          #                           "2.5%" = vy, 
          #                           "50%" = vy, 
          #                           "97.5%" = vy))
          break
        } else {
          rxxw_ML <- rbind(rxxw_ML, 
                           data.frame(temperatures=temperature, 
                                      Mean = vy, sd = NA, 
                                      "2.5%" = ic5[1, 2], 
                                      "50%" = vy, 
                                      "97.5%" = ic5[1, 3]))
        }
      }
      if (warn) {
        # J'ai une erreur avec la méthode delta
        
        sigma <- solve(hessian)
        s. <- svd(sigma)
        R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
        
        pre0.9_9994 = FALSE
        
        dataR <- matrix(rnorm(replicate.CI * ncol(sigma)), nrow = replicate.CI, byrow = !pre0.9_9994) %*% R
        dataR <- sweep(dataR, 2, parameters[rownames(hessian)], "+")
        colnames(dataR) <- rownames(hessian)
        
        ajouter <- matrix(rep(fixed.parameters, replicate.CI), 
                          nrow=replicate.CI, byrow=TRUE,
                          dimnames = list(c(NULL), names(fixed.parameters)))
        dataR <- cbind(dataR, ajouter)
        
        xxw <- matrix(ncol = 2, 
                      nrow = length(temp)*nrow(dataR))
        xxw2 <- matrix(data = as.numeric(), ncol = length(temp), 
                       nrow = nrow(dataR))
        
        xxw[, 1] <- rep(temp, nrow(dataR))
        for (i in 1:nrow(dataR)) {
          r <- SSM(T=temp, 
                   parms=c(dataR[i,], fixed.parameters))[[set.par]]
          jhn <- (1+(i-1)*length(temp)):((i-1)*length(temp)+length(temp))
          # xxw[jhn, 1] <- temp
          xxw[jhn, 2] <- r
          xxw2[i, ] <- r
        }
        rxxw_ML <- apply(xxw2, MARGIN = 2, 
                      FUN = function(x) c(Mean=mean(x), sd=sd(x), quantile(x, probs =c((1-probs)/2, 0.5, 1-(1-probs)/2), na.rm=TRUE)))
        rxxw_ML <- t(rbind(temperatures=temp, rxxw_ML))
      }
    rxxw_ML <- as.matrix(t(rxxw_ML))
    }
    
      
    # Si je n'ai pas de tableau generes, je les mets a NA
    if (is.null(rxxw))
          rxxw <- as.matrix(t(data.frame(temperatures=temp, Mean=NA, sd=NA, "2.5%"=NA, "50%"=NA, "97.5%"=NA)))
    if (is.null(rxxw_ML))
          rxxw_ML <- as.matrix(t(data.frame(temperatures=temp, Mean=NA, sd=NA, "2.5%"=NA, "50%"=NA, "97.5%"=NA)))
    
    
    if (!is.null(c(parameters, fixed.parameters))) {
      r <- SSM(T=temp, 
               parms=c(parameters, fixed.parameters))[[set.par]]
      r[is.infinite(r)] <- NA
      xxw_ML <- data.frame(temperatures=temp, r=r)
      rxxw_ML <- rbind(rxxw_ML, ML=r)
    }  else {
      xxw_ML <- data.frame(temperatures=temp, r=NA)
      rxxw_ML <- rbind(rxxw_ML, ML=NA)
    }
    
    if (is.null(xxw)) xxw <- xxw_ML
    if (scaleY=="auto") scaleY <- 10^(-floor(log10(max(xxw[, 2], na.rm = TRUE))))
    
    xxw[, 2] <- xxw[, 2] * scaleY
    
    rxxw[2:6, ] <- rxxw[2:6, ] * scaleY
    rxxw_ML[2:7, ] <- rxxw_ML[2:7, ] * scaleY
    rownames(rxxw) <- c("temperatures", "Mean", "sd", "X2.5", "X50", "X97.5")
    rownames(rxxw_ML) <- c("temperatures", "Mean", "sd", "X2.5", "X50", "X97.5", "ML")
    
    if (is.null(ylab)) ylab <- as.expression(bquote(.(paste0("r x ", format(scaleY, scientific=FALSE), 
                                                             " (mm.min"))*""^-1*")"))
    
    if (show.density & !is.null(xxw)) {
      if (new) {
        par(new=FALSE)
    smoothScatter(x=xxw[, 1], y=xxw[, 2], las=las, bty = bty,  pch=pch, 
                  bandwidth = bandwidth, 
                  xlab=xlab, 
                  main=main,
                  colramp=colramp, 
                  ylab=ylab, 
                  ylim=ylim, xlim = xlim, nbin = 128, postPlotHook=NULL, 
                  cex.axis=cex.axis, 
                  cex.lab=cex.lab)
      } else {
        par(new=TRUE)
        smoothScatter(x=xxw[, 1], y=xxw[, 2], bty = "n",  pch=pch, 
                      bandwidth = bandwidth, 
                      xlab="", 
                      main="",
                      colramp=colramp, 
                      ylab="", 
                      # xlim=c(ScalePreviousPlot()$xlim[1]-ScalePreviousPlot()$xlim[4]*0.04, ScalePreviousPlot()$xlim[2]+ScalePreviousPlot()$xlim[4]*0.04), 
                      # ylim=c(ScalePreviousPlot()$ylim[1]-ScalePreviousPlot()$ylim[4]*0.04, ScalePreviousPlot()$ylim[2]+ScalePreviousPlot()$ylim[4]*0.04),
                      ylim = ScalePreviousPlot()$ylim[1:2], 
                      xlim = ScalePreviousPlot()$xlim[1:2], 
                      nbin = 128, postPlotHook=NULL, 
                      axes=FALSE,
                      new=FALSE)
      }
    } else {
      if (new) {
        par(new=FALSE)
      plot(xxw[, 1], y=xxw[, 2], type="n", 
           las=las, 
           bty = bty, 
           xlab=xlab, 
           main=main,
           ylab=ylab, 
           ylim=ylim, 
           xlim = xlim, 
           cex.axis=cex.axis, 
           cex.lab=cex.lab)
      # } else {
      #   par(new=TRUE)
      #   plot(xxw[, 1], y=xxw[, 2], type="n", 
      #        las=las, 
      #        bty = bty, 
      #        xlab="", 
      #        main="",
      #        ylab="", 
      #        ylim = ScalePreviousPlot()$ylim[1:2], 
      #        xlim = ScalePreviousPlot()$xlim[1:2], 
      #        axes=FALSE)
      }
    }
    

    vx <- c(rxxw["temperatures", ], rev(rxxw["temperatures", ]))
    
    if (any(curves == "mcmc quantiles") & any(!is.na(rxxw[c("X2.5", "X50", "X97.5"), ]))) {
      if (polygon) {
        vy <- c(rxxw["X2.5", ], rev(rxxw["X97.5", ]))
        vy <- ifelse(vy<ylim[1], ylim[1], vy)
        vy <- ifelse(vy>ylim[2], ylim[2], vy)
        polygon(x=vx, y=vy, col=col.polygon, border = NA)
      }
      
    dxy <- data.frame(x=rxxw["temperatures", ], y=ifelse(rxxw["X50", ]>ylim[1] & rxxw["X50", ]<ylim[2], rxxw["X50", ], NA))
    dxy <- na.omit(dxy)
    lines(x = dxy$x, y=dxy$y, lty=lty, lwd=lwd, col=col)
    
    dxy <- data.frame(x=rxxw["temperatures", ], y=ifelse(rxxw["X2.5", ]>ylim[1] & rxxw["X2.5", ]<ylim[2], rxxw["X2.5", ], NA))
    dxy <- na.omit(dxy)
    lines(x = dxy$x, y=dxy$y, lty=ltyCI, lwd=lwdCI, col=col)
    
    dxy <- data.frame(x=rxxw["temperatures", ], y=ifelse(rxxw["X97.5", ]>ylim[1] & rxxw["X97.5", ]<ylim[2], rxxw["X97.5", ], NA))
    dxy <- na.omit(dxy)
    lines(x = dxy$x, y=dxy$y, lty=ltyCI, lwd=lwdCI, col=col)
    }
    
    if (any(curves == "mcmc mean-sd") & any(!is.na(rxxw[c("Mean", "sd"), ]))) {
      if (polygon) {
        vy <- c(rxxw["Mean", ]+1.96*rxxw["sd", ], rev(rxxw["Mean", ]-1.96*rxxw["sd", ]))
        vy <- ifelse(vy<ylim[1], ylim[1], vy)
        vy <- ifelse(vy>ylim[2], ylim[2], vy)
        polygon(x=vx, y=vy, col=col.polygon, border = NA)
      }
      
      dxy <- data.frame(x=rxxw["temperatures", ], y=ifelse(rxxw["Mean", ]>ylim[1] & rxxw["Mean", ]<ylim[2], rxxw["Mean", ], NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=lty, lwd=lwd, col=col)
      
      dxy <- data.frame(x=rxxw["temperatures", ], y=ifelse((rxxw["Mean", ]-1.96*rxxw["sd", ])>ylim[1] & (rxxw["Mean", ]-1.96*rxxw["sd", ])<ylim[2], rxxw["Mean", ]-1.96*rxxw["sd", ], NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=ltyCI, lwd=lwdCI, col=col)
      
      dxy <- data.frame(x=rxxw["temperatures", ], y=ifelse((rxxw["Mean", ]+1.96*rxxw["sd", ])>ylim[1] & (rxxw["Mean", ]+1.96*rxxw["sd", ])<ylim[2], rxxw["Mean", ]+1.96*rxxw["sd", ], NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=ltyCI, lwd=lwdCI, col=col)
    }
    if (any(curves == "ml") & !is.null(parameters) & any(!is.na(rxxw_ML[c("ML"), ]))) {
      dxy <- data.frame(x=rxxw_ML["temperatures", ], y=ifelse(rxxw_ML["ML", ]>ylim[1] & rxxw_ML["ML", ]<ylim[2], rxxw_ML["ML", ], NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=lty, lwd=lwd, col=col)
    }
    
    if (any(curves == "ml quantiles") & !is.null(parameters) & any(!is.na(rxxw_ML[c("X2.5", "X50", "X97.5"), ]))) {
      if (polygon) {
        vy <- c(rxxw_ML["X2.5", ], rev(rxxw_ML["X97.5", ]))
        vy <- ifelse(vy<ylim[1], ylim[1], vy)
        vy <- ifelse(vy>ylim[2], ylim[2], vy)
        polygon(x=vx, y=vy, col=col.polygon, border = NA)
      }
      
      dxy <- data.frame(x=rxxw_ML["temperatures", ], y=ifelse(rxxw_ML["X50", ]>ylim[1] & rxxw_ML["X50", ]<ylim[2], rxxw_ML["X50", ], NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=lty, lwd=lwd, col=col)
      
      dxy <- data.frame(x=rxxw_ML["temperatures", ], y=ifelse(rxxw_ML["X2.5", ]>ylim[1] & rxxw_ML["X2.5", ]<ylim[2], rxxw_ML["X2.5", ], NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=ltyCI, lwd=lwdCI, col=col)
      
      dxy <- data.frame(x=rxxw_ML["temperatures", ], y=ifelse(rxxw_ML["X97.5", ]>ylim[1] & rxxw_ML["X97.5", ]<ylim[2], rxxw_ML["X97.5", ], NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=ltyCI, lwd=lwdCI, col=col)
    }
    
    if (!is.null(temperatures) & show.hist) {
    
    xlim <- ScalePreviousPlot()$xlim[c("begin", "end")]
    ylim <- ScalePreviousPlot()$ylim[c("begin", "end")]
    
    par(new=TRUE)
    if (log.hist==FALSE) {
     if (is.null(ylimH)) {
       a <- hist(x=temperatures, xlab="", ylab="", main="", axes=FALSE, freq=FALSE, 
                            xlim=xlim, breaks=breaks)
      } else {
        a <- hist(x=temperatures, xlab="", ylab="", main="", axes=FALSE, freq=FALSE, 
                  xlim=xlim, ylim=ylimH, breaks=breaks)
      }
      axis(side=4, las=1, at=atH, cex.axis=cex.axis)
      
    } else {
      ax <- hist(x=temperatures, plot=FALSE)
      yi <- ifelse(log(ax$histogram$density*1000)<0, 0, log(ax$histogram$density*1000))
      myi <- max(yi)
      yi <- yi/myi
      if (!is.null(ylimH)) yi <- yi*ylimH[2]
      
      yi <- yi*ylim[2]
      
      for (i in seq_along(ax$histogram$density)) {
        polygon(x=ax$histogram$breaks[c(i, i+1, i+1, i)], y=c(0, 0, yi[i], yi[i]))
      }
      
      if (is.null(atH)) {
        atH <- seq(from=0, to=max(ax$histogram$density), by=0.05)
      } 
        
        yi <- ifelse(log(atH*1000)<0, 0, log(atH*1000))
        # myi <- max(yi)
        yi <- yi/myi
        
        if (!is.null(ylimH)) yi <- yi*ylimH[2]
        yi <- yi*ylim[2]
        
        axis(side=4, las=1, at=yi, labels = atH, cex.axis=cex.axis)
      
        atH <- yi
    }
    
    mtext(ylabH, side=4, line=3, at=ifelse(is.null(atH), NA, mean(c(atH[1], rev(atH)[1]))), cex=cex.lab)
    
    par(new=TRUE)
    # je retablis l'echelle des y et celle de R
    plot(x = 1, y=1, ylim=ylim, xlim=xlim, xlab="", ylab="", axes=FALSE, bty="n", type="n")
    
    }
    
    return(invisible(scaleY))
  }
