#' plotR plots the fitted growth rate dependent on temperature and the density of the mcmc
#' @title Plot the fitted growth rate dependent on temperature and its density
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with the value of scaleY to be used with other plotR function and the plot data in xy list element
#' @param result A result object or a list of result objects
#' @param resultmcmc A result object from GRTN_MHmcmc() function
#' @param chain The chain to use in resultmcmc
#' @param parameters A set of parameters - Has the priority over result
#' @param fixed.parameters A set of fixed parameters
#' @param temperatures A set of temperatures - Has the priority over result
#' @param hessian An hessian matrix
#' @param replicate.CI Number of replicates to estimate confidence interval with Hessian if delta method failed
#' @param curve What curve to show: "MCMC quantiles" or "MCMC mean-SD" based on mcmc or "ML" or "ML quantiles" or "ML mean-SE" for maximum-likelihood. Or "none"
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
#' @param yaxt The yaxt parameter of y-axis
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
#' The parameter curve is case insensitive. If only parameters is given, curve must be ML.
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' # Note that the confidence interval is not the same for mcmc and ml quantiles
#' plotR(result = resultNest_4p_SSM, 
#'              resultmcmc=resultNest_mcmc_4p_SSM, 
#'              curve = "MCMC quantiles", ylim=c(0, 8))
#' plotR(resultNest_4p_SSM, curve="ML quantiles", ylim=c(0, 6))
#' #################
#' plotR(resultmcmc=resultNest_mcmc_4p_SSM, ylim=c(0, 10), 
#'              curve = "MCMC quantiles", show.density=TRUE)
#' #################
#' plotR(resultmcmc=resultNest_mcmc_4p_SSM, 
#'              curve = "MCMC quantiles", polygon=TRUE)
#' #################
#' plotR(resultmcmc=resultNest_mcmc_6p_SSM, ylim=c(0,8), 
#'       curve = "MCMC quantiles", polygon=TRUE, col.polygon = rgb(0, 1, 0, 1))
#' plotR(resultmcmc=resultNest_mcmc_4p_SSM, ylim=c(0,8), 
#'        curve = "MCMC quantiles", polygon=TRUE, col.polygon = rgb(1, 0, 0, 0.5), 
#'        new=FALSE)
#' legend("topleft", legend=c("SSM 4 parameters", "SSM 6 parameters"), 
#'         pch=c(15, 15), col=c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 1)))
#' #################
#' sy <- plotR(resultmcmc=resultNest_mcmc_4p_SSM, ylim=c(0, 8), 
#'              curve = "MCMC quantiles", show.density=FALSE)
#' plotR(resultmcmc=resultNest_mcmc_6p_SSM, col="red", ylim=c(0, 8), 
#'              curve = "MCMC quantiles", show.density=FALSE, 
#'              new=FALSE, scaleY=sy)
#' #################
#' sy <- plotR(result=resultNest_6p_SSM, curve="none", 
#'              scaleY=1E5, 
#'              ylim=c(0, 8), 
#'              show.hist = TRUE, new = TRUE, mar=c(4, 4, 1, 4))
#' #################
#' plotR(result=resultNest_6p_SSM, curve="ML", 
#'              ylim=c(0, 8), 
#'              show.hist = TRUE, ylimH=c(0,1), atH=c(0, 0.1, 0.2))
#' ################
#' plotR(result = resultNest_4p_SSM, ylim=c(0, 8), 
#'              resultmcmc=resultNest_mcmc_4p_SSM, 
#'              show.density = TRUE, 
#'              curve = "MCMC quantiles")
#' #################
#' plotR(resultmcmc=resultNest_mcmc_4p_SSM, ylim=c(0, 8), 
#'              curve = "MCMC quantiles", show.density=TRUE, scaleY=1E5)
#' }
#' @export


plotR <-
  function(result = NULL                                                                   , 
           resultmcmc = NULL                                                               , 
           chain=1                                                                         , 
           parameters = NULL                                                               , 
           fixed.parameters = NULL                                                         , 
           temperatures  = NULL                                                            ,
           curve = "ML quantiles"                                                          , 
           set.par=1                                                                       , 
           ylim=c(0, 5)                                                                    , 
           xlim=c(20,35)                                                                   , 
           hessian = NULL                                                                  , 
           replicate.CI = 1000                                                             , 
           cex.lab = par("cex")                                                            , 
           cex.axis = par("cex")                                                           ,
           scaleY="auto"                                                                   , 
           lty=1                                                                           , 
           ltyCI=3                                                                         , 
           lwd=1                                                                           , 
           lwdCI=1                                                                         , 
           col = "black"                                                                   , 
           col.polygon="grey"                                                              , 
           polygon=FALSE                                                                   , 
           probs=0.95                                                                      ,
           colramp=colorRampPalette(c("white", rgb(red = 0.5, green = 0.5, blue = 0.5)))   , 
           bandwidth = c(0.1, 0.01)                                                        , 
           pch = ""                                                                        , 
           main=""                                                                         ,
           xlab = expression("Temperature in "*degree*"C")                                 , 
           ylab = NULL                                                                     , 
           yaxt = "s"                                                                       ,
           bty = "n"                                                                       , 
           las = 1                                                                         , 
           by.temperature=0.1                                                              , 
           show.density=FALSE                                                              , 
           new=TRUE                                                                        , 
           show.hist=FALSE                                                                 , 
           ylimH = NULL                                                                    , 
           atH = NULL                                                                      , 
           ylabH="Temperature density"                                                     , 
           breaks = "Sturges"                                                              , 
           log.hist=FALSE                                                                  ,
           mar=NULL                                                                        ) {
    
    # result = NULL; resultmcmc=NULL; chain=1; parameters = NULL; fixed.parameters = NULL; hessian = NULL; replicate.CI=1000; probs=0.95; temperatures  = NULL;curve = "ML quantiles"; set.par=1; ylim=c(0, 5); xlim=c(20,35); cex.lab = 1; cex.axis = 1;scaleY="auto"; lty=1; ltyCI=3; lwd=1; lwdCI=1; colramp=colorRampPalette(c("white", rgb(red = 0.5, green = 0.5, blue = 0.5))); bandwidth = c(0.3, 0.05); pch = ""; main=""; col = "black"; col.polygon="grey"; polygon=FALSE; xlab = expression("Temperature in"*degree*"C"); ylab = NULL; bty = "n"; las = 1; by.temperature=0.1; show.density=FALSE;new=TRUE;show.hist=FALSE; ylimH = NULL; atH = NULL;ylabH="Temperature density";breaks = "Sturges";log.hist=FALSE;mar=NULL
    # resultmcmc = resultNest_mcmc_4p_SSM
    # result = resultNest_4p_SSM
    
    return_Out <- NULL
    
    curve <- tolower(curve)
    curve <- match.arg(curve, choices = c("mcmc quantiles", "mcmc mean-sd", "ml", "ml quantiles", "ml mean-se", "none"))
    
    SSM <- getFromNamespace(".wrapperSSM", ns="embryogrowth")
    
    if (inherits(result, "mcmcComposite")) {
      resultmcmc <- result
      result <- NULL
    }
    
    if (inherits(result, "numeric")) {
      parameters <- result
      result <- NULL
    }
    
    if (is.null(hessian) & !is.null(result)) hessian <- result$hessian
    if (is.null(fixed.parameters) & !is.null(result)) fixed.parameters <- result$fixed.parameters
    if (is.null(parameters) & !is.null(result)) parameters <- result$par
    
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
      if (!is.null(temperatures) & show.hist & new) {
        par(mar=c(4, 4, 1, 4)+0.4) 
      } else {
        par(mar=c(4, 4, 4, 1)+0.4)
      }
    } else {
      if (all(!is.na(mar))) par(mar=mar)
    }
    
    temp <- seq(from=xlim[1], to=xlim[2], by=by.temperature)
    RandomMatrixRforTempsAndReplicates <- NULL
    
    if (((curve=="ml quantiles") | (curve == "ml mean-se")) & is.null(hessian)) {
      curve <- "ml"
      replicate.CI <- 0
      warning("ml quantiles or mean-se curve needs an hessian matrix")
    }
    
    if (!is.null(resultmcmc) & (show.density | curve=="mcmc quantiles" | curve=="mcmc mean-sd")) {
      RandomMatrixRforTempsAndReplicates <- RandomFromHessianOrMCMC(mcmc = resultmcmc, replicates = replicate.CI, 
                                                                    method="mcmc", 
                                                                    fixed.parameters = fixed.parameters, 
                                                                    chain = chain, 
                                                                    ParTofn = "parms", 
                                                                    fn=SSM, 
                                                                    set.par=set.par, 
                                                                    T=temp)
    }
    
    # j'ai un hessian qui est fourni
    if (curve=="ml quantiles" & !is.null(hessian)) {
      RandomMatrixRforTempsAndReplicates <- RandomFromHessianOrMCMC(Hessian = hessian                       , 
                                                                    replicates = replicate.CI               , 
                                                                    fitted.parameters = parameters          , 
                                                                    method = "Hessian"                      , 
                                                                    fixed.parameters = fixed.parameters     , 
                                                                    ParTofn = "parms"                       , 
                                                                    fn=SSM                                  , 
                                                                    set.par=set.par                         ,  
                                                                    T=temp                                  )
    }
    
    if (curve=="ml mean-se" & !is.null(hessian)) {
      RandomMatrixRforTempsAndReplicates <- RandomFromHessianOrMCMC(se = SEfromHessian(hessian)        , 
                                                                    replicates = replicate.CI          , 
                                                                    fitted.parameters = parameters     , 
                                                                    method = "SE"                      , 
                                                                    fixed.parameters = fixed.parameters, 
                                                                    ParTofn = "parms"                  , 
                                                                    fn=SSM                             , 
                                                                    set.par=set.par                    ,  
                                                                    T=temp                             )
    }
    
    if (is.null(RandomMatrixRforTempsAndReplicates))
      if (!is.null(c(parameters, fixed.parameters)) & curve=="ml") {
        RandomMatrixRforTempsAndReplicates <- RandomFromHessianOrMCMC(replicates = 0                     , 
                                                                      fitted.parameters = parameters     , 
                                                                      fixed.parameters = fixed.parameters, 
                                                                      method="null"                      , 
                                                                      ParTofn = "parms"                  , 
                                                                      fn=SSM                             , 
                                                                      set.par=set.par                    , 
                                                                      T=temp                             )
        
      }
    
    if (is.null(RandomMatrixRforTempsAndReplicates)) {
      MatrixRforSmooth <- NULL
      MatrixRforTempsAndReplicates <- NULL
      SynthesisMatrixRforTempsAndReplicates <- NULL
    } else {
      
      # MatrixRforTempsAndReplicates <- RandomMatrixRforTempsAndReplicates$fn[, (1:length(temp))+(set.par-1)*length(temp), drop=FALSE]
      # 19-04-2021: set par est maintenant inclus dans RandomMatrixRforTempsAndReplicates
      MatrixRforTempsAndReplicates <- RandomMatrixRforTempsAndReplicates$fn
      
      # Je crée une matrice avec les moyennes, sd et quantiles pour chaque temp
      SynthesisMatrixRforTempsAndReplicates <- apply(MatrixRforTempsAndReplicates, MARGIN = 2, 
                                                     FUN = function(x) c(Mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE), quantile(x, probs =c((1-probs)/2, 0.5, 1-(1-probs)/2), na.rm=TRUE)))
      # Et je rajoute une ligne avec les temp
      SynthesisMatrixRforTempsAndReplicates <- rbind(temperatures=temp, SynthesisMatrixRforTempsAndReplicates)
      
      
      MatrixRforSmooth <- matrix(data=as.numeric(), ncol = 2, 
                                 nrow = length(temp)*replicate.CI)
      MatrixRforSmooth[, 1] <- rep(temp, replicate.CI)
      MatrixRforSmooth[, 2] <- as.vector(t(MatrixRforTempsAndReplicates))
      MatrixRforSmooth <- na.omit(MatrixRforSmooth)
    }
    
    if ((scaleY=="auto") & (is.null(MatrixRforSmooth))) {
      stop("scaleY cannot be estimated")
    }
    
    if (scaleY=="auto") scaleY <- 10^(-floor(log10(max(MatrixRforSmooth[is.finite(MatrixRforSmooth[, 2]), 2], na.rm = TRUE))))
    
    MatrixRforSmooth[, 2] <- MatrixRforSmooth[, 2] * scaleY
    
    SynthesisMatrixRforTempsAndReplicates[2:6, ] <- SynthesisMatrixRforTempsAndReplicates[2:6, ] * scaleY
    if (!is.null(SynthesisMatrixRforTempsAndReplicates)) {
      rownames(SynthesisMatrixRforTempsAndReplicates) <- c("temperatures", "Mean", "sd", "X2.5", "X50", "X97.5")
    }
    if (is.null(ylab)) ylab <- as.expression(bquote(.(paste0("r x ", format(scaleY, scientific=FALSE), 
                                                             " (min"))*""^-1*")"))
    
    if (show.density & !is.null(MatrixRforSmooth)) {
      if (new) {
        par(new=FALSE)
        smoothScatter(x=MatrixRforSmooth[, 1], y=MatrixRforSmooth[, 2], las=las, bty = bty,  pch=pch, 
                      bandwidth = bandwidth, 
                      xlab=xlab, 
                      main=main,
                      colramp=colramp, 
                      ylab=ylab, 
                      ylim=ylim, xlim = xlim, nbin = 128, postPlotHook=NULL, 
                      cex.axis=cex.axis, 
                      cex.lab=cex.lab, 
                      yaxt=yaxt)
      } else {
        par(new=TRUE)
        smoothScatter(x=MatrixRforSmooth[, 1], y=MatrixRforSmooth[, 2], bty = "n",  pch=pch, 
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
                      new=FALSE, 
                      yaxt=yaxt)
      }
    } else {
      if (new) {
        par(new=FALSE)
        plot(MatrixRforSmooth[, 1], y=MatrixRforSmooth[, 2], type="n", 
             las=las, 
             bty = bty, 
             xlab=xlab, 
             main=main,
             ylab=ylab, 
             ylim=ylim, 
             xlim = xlim, 
             cex.axis=cex.axis, 
             cex.lab=cex.lab, 
             yaxt=yaxt)
      }
    }
    
    
    if (!is.null(SynthesisMatrixRforTempsAndReplicates)) {
      vx <- c(SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
              rev(SynthesisMatrixRforTempsAndReplicates["temperatures", ]))
    } else {
      vx <- c(temp, rev(temp))
    }
    
    if ((curve == "ml quantiles") | (curve == ("mcmc quantiles"))) {
      return_Out <- data.frame(x=SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
                               y=SynthesisMatrixRforTempsAndReplicates["X50", ], 
                               y.lower=SynthesisMatrixRforTempsAndReplicates["X2.5", ], 
                               y.upper=SynthesisMatrixRforTempsAndReplicates["X97.5", ])
      if (polygon) {
        vy <- c(SynthesisMatrixRforTempsAndReplicates["X2.5", ], rev(SynthesisMatrixRforTempsAndReplicates["X97.5", ]))
        vy <- ifelse(vy<ylim[1], ylim[1], vy)
        vy <- ifelse(vy>ylim[2], ylim[2], vy)
        polygon(x=vx, y=vy, col=col.polygon, border = NA)
      }
      
      dxy <- data.frame(x=SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
                        y=ifelse(SynthesisMatrixRforTempsAndReplicates["X50", ]>=ylim[1] & SynthesisMatrixRforTempsAndReplicates["X50", ]<=ylim[2], 
                                 SynthesisMatrixRforTempsAndReplicates["X50", ], 
                                 NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=lty, lwd=lwd, col=col)
      
      dxy <- data.frame(x=SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
                        y=ifelse(SynthesisMatrixRforTempsAndReplicates["X2.5", ]>=ylim[1] & SynthesisMatrixRforTempsAndReplicates["X2.5", ]<=ylim[2], 
                                 SynthesisMatrixRforTempsAndReplicates["X2.5", ], 
                                 NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=ltyCI, lwd=lwdCI, col=col)
      
      dxy <- data.frame(x=SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
                        y=ifelse(SynthesisMatrixRforTempsAndReplicates["X97.5", ]>=ylim[1] & SynthesisMatrixRforTempsAndReplicates["X97.5", ]<=ylim[2], 
                                 SynthesisMatrixRforTempsAndReplicates["X97.5", ], 
                                 NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=ltyCI, lwd=lwdCI, col=col)
    }
    
    if ((curve == "mcmc mean-sd") | (curve == "ml mean-se")) {
      return_Out <- data.frame(x=SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
                               y=SynthesisMatrixRforTempsAndReplicates["Mean", ], 
                               y.lower=SynthesisMatrixRforTempsAndReplicates["Mean", ]-1.96*SynthesisMatrixRforTempsAndReplicates["sd", ], 
                               y.upper=SynthesisMatrixRforTempsAndReplicates["Mean", ]+1.96*SynthesisMatrixRforTempsAndReplicates["sd", ])
      
      if (polygon) {
        vy <- c(SynthesisMatrixRforTempsAndReplicates["Mean", ]+1.96*SynthesisMatrixRforTempsAndReplicates["sd", ], rev(SynthesisMatrixRforTempsAndReplicates["Mean", ]-1.96*SynthesisMatrixRforTempsAndReplicates["sd", ]))
        vy <- ifelse(vy<ylim[1], ylim[1], vy)
        vy <- ifelse(vy>ylim[2], ylim[2], vy)
        polygon(x=vx, y=vy, col=col.polygon, border = NA)
      }
      
      dxy <- data.frame(x=SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
                        y=ifelse(SynthesisMatrixRforTempsAndReplicates["Mean", ]>=ylim[1] & SynthesisMatrixRforTempsAndReplicates["Mean", ]<=ylim[2], 
                                 SynthesisMatrixRforTempsAndReplicates["Mean", ], 
                                 NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=lty, lwd=lwd, col=col)
      
      dxy <- data.frame(x=SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
                        y=ifelse((SynthesisMatrixRforTempsAndReplicates["Mean", ]-1.96*SynthesisMatrixRforTempsAndReplicates["sd", ])>=ylim[1] & (SynthesisMatrixRforTempsAndReplicates["Mean", ]-1.96*SynthesisMatrixRforTempsAndReplicates["sd", ])<=ylim[2], SynthesisMatrixRforTempsAndReplicates["Mean", ]-1.96*SynthesisMatrixRforTempsAndReplicates["sd", ], NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=ltyCI, lwd=lwdCI, col=col)
      
      dxy <- data.frame(x=SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
                        y=ifelse((SynthesisMatrixRforTempsAndReplicates["Mean", ]+1.96*SynthesisMatrixRforTempsAndReplicates["sd", ])>=ylim[1] & (SynthesisMatrixRforTempsAndReplicates["Mean", ]+1.96*SynthesisMatrixRforTempsAndReplicates["sd", ])<=ylim[2], SynthesisMatrixRforTempsAndReplicates["Mean", ]+1.96*SynthesisMatrixRforTempsAndReplicates["sd", ], NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=ltyCI, lwd=lwdCI, col=col)
    }
    
    if (curve == "ml") {
      return_Out <- data.frame(x=SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
                               y=SynthesisMatrixRforTempsAndReplicates["Mean", ], 
                               y.lower=NA, 
                               y.upper=NA)
      
      
      dxy <- data.frame(x=SynthesisMatrixRforTempsAndReplicates["temperatures", ], 
                        y=ifelse(SynthesisMatrixRforTempsAndReplicates["Mean", ]>=ylim[1] & SynthesisMatrixRforTempsAndReplicates["Mean", ]<=ylim[2], SynthesisMatrixRforTempsAndReplicates["Mean", ], NA))
      dxy <- na.omit(dxy)
      lines(x = dxy$x, y=dxy$y, lty=lty, lwd=lwd, col=col)
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
      plot(x = 1, y=1, ylim=ylim, xlim=xlim, xlab="", ylab="", axes=FALSE, bty="n", type="n", yaxt=yaxt)

    }
    
    return(invisible(list(scaleY=scaleY, xy=return_Out)))
  }
