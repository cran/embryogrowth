#' plotR shows the fitted growth rate dependent on temperature
#' @title Show the fitted growth rate dependent on temperature
#' @author Marc Girondot
#' @return A list with data.frame with the confidence interval and the average.
#' @param result A result object or a list of result objects
#' @param ... Parameters for plot() such as main= or ylim=
#' @param parameters Indicate some parameters if the result object is not supplied
#' @param fixed.parameters Indicate some parameters if the result object is not supplied
#' @param SE The standard error for the parameters or a list of SE if several results. Use NA to force not use SE
#' @param x.SE The factor used to show the confidence interval envelope
#' @param set.par 1 or 2 or a list of 1 or 2 to designate with set of parameters to show
#' @param size If indicated, will show the growth rate for this size. Useful only for model with two sets of parameters, High and Low and a transition
#' @param legend Text to show in bottom right legend or a list of text if several results
#' @param show.legend Should the legend about several series be shown?
#' @param col The color to use for a list of colors if several results
#' @param lty The type of line to use if several results as a list
#' @param ltyCI The type of line to use for confidence interval as a list
#' @param lwd The type of line to use if several results as a list
#' @param lwdCI The type of line to use for confidence interval as a list
#' @param xlim Range of values for x-axis
#' @param xlimR Range of values to be displayed for R curve or vector of values; can be a list if a list of results is used
#' @param xlimSE Range of values to be displayed for SE curves; can be a list if a list of results is used
#' @param scaleY Scaling factor for y axis or "auto"
#' @param replicate.CI Number of replicates to estimate CI
#' @param show.box If TRUE show a box with "mean" and "confidence interval"
#' @param local.box Position of the box with "mean" and "confidence interval", default="topleft"
#' @param pch.anchors Symbol used to show anchors
#' @param cex.anchors Size of symbol used to show anchors
#' @param col.anchors Color of symbols used to show anchors
#' @param y.anchors Position of anchors in y axis
#' @param show.anchors Should the anchors been shown
#' @description To show the growth rate, the syntaxe is:\cr
#' plotR(result=res)\cr
#' If SE is a matrix with two raws, the row values will use to draw the CI.
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' resultNest_4p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p)
#' plotR(resultNest_4p, ylim=c(0,3))
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p, accept=TRUE)
#' # Take care, it can be very long, sometimes several days
#' result_mcmc_4p_80 <- GRTRN_MHmcmc(result=resultNest_4p,  
#' 	parametersMCMC=pMCMC, n.iter=10000, n.chains = 1, n.adapt = 0,  
#' 	thin=1, trace=TRUE)
#' data(result_mcmc_4p)
#' plotR(result=resultNest_4p, SE=result_mcmc_4p$SD,  
#'  ylim=c(0, 3), x.SE=1)
#' x <- structure(c(115.758929130522, 428.649022170996, 503.687251738993, 
#' 12.2621455821612, 306.308841227278, 116.35048615105), .Names = c("DHA", 
#' "DHH", "DHL", "DT", "T12L", "Rho25"))
#' plotR(parameters=x, xlim=c(20, 35))
#' pfixed <- c(rK=2.093313)
#' resultNest_6p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_6p)
#' plotR(list(resultNest_4p, resultNest_6p), ylim=c(0, 3), 
#' col=c("black", "red"), legend=c("4 parameters", "6 parameters"))
#' ##########################################
#' # new formulation of parameters using anchors
#' data(resultNest_newp)
#' # without envelope
#' plotR(resultNest_newp, ylim=c(0, 5))
#' # with envelope based in 1.96*SE and central curve based on mean
#' plotR(result=resultNest_newp, ylim=c(0, 5), 
#'  SE=result_mcmc_newp$SD)
#' # with envelope based on quantiles and central curve based on mean
#' plotR(result=resultNest_newp, ylim=c(0, 5), 
#'  SE=apply(result_mcmc_newp[["resultMCMC"]][[1]], 
#'    MARGIN=2, FUN=quantile, probs=c(0.025, 0.975)))
#' # with envelope based on quantiles and central curve based on median
#' plotR(result=resultNest_newp, ylim=c(0, 5), 
#'  SE=apply(result_mcmc_newp[["resultMCMC"]][[1]], 
#'    MARGIN=2, FUN=quantile, probs=c(0.025, 0.975)),
#'  parameters=apply(result_mcmc_newp[["resultMCMC"]][[1]], 
#'    MARGIN=2, FUN=quantile, probs=c(0.5)))
#' # Example to get the results
#' (plotR(result=resultNest_newp, ylim=c(0, 5), 
#'        SE=apply(result_mcmc_newp[["resultMCMC"]][[1]], 
#'                 MARGIN=2, FUN=quantile, probs=c(0.025, 0.975)),
#'        parameters=apply(result_mcmc_newp[["resultMCMC"]][[1]],
#'                        MARGIN=2, FUN=quantile, probs=c(0.5)), 
#'        xlimR=as.numeric(names(resultNest_newp$par))-273.15)[[1]])
#'  ##########################################
#'  # New Weilbull model
#'  ##########################################
#'  x <- c(k=3, lambda=3, theta=290)
#'  resultNest_4p_Weibull <- searchR(parameters=x, fixed.parameters=pfixed,
#'    temperatures=formated, derivate=dydt.Gompertz, M0=1.7,
#'    test=c(Mean=39.33, SD=1.92))
#'  plotR(resultNest_4p_Weibull)

#' }
#' @export


plotR <-
  function(result=NULL, parameters=NULL, fixed.parameters=NULL, col="black", legend=NA, 
           show.legend=TRUE,
           SE=NULL, x.SE=qnorm(0.975), set.par=1, size=NULL, xlim=c(20,35), 
           scaleY="auto", lty=1, ltyCI=3, lwd=1, lwdCI=1, 
           xlimR=xlim, xlimSE=xlim, replicate.CI=100, show.box=TRUE, local.box="topleft", 
           pch.anchors=19, cex.anchors=1, col.anchors="black", y.anchors=0, 
           show.anchors=TRUE, ...) {
    
    # result=NULL;parameters=NULL; fixed.parameters=NULL; lty=1; ltyCI=3; lwd=1; lwdCI=1; col="black"; legend=NA; SE=NULL; set.par=1; size=NA; xlim=c(20,35); xlimR=xlim; xlimSE=xlim; x.SE=qnorm(0.975); scaleY="auto"; replicate.CI=100; show.box=TRUE; local.box="topleft"; pch.anchors=19; cex.anchors=1; col.anchors="black"; y.anchors=0; show.anchors=TRUE
    # result <- resultNest_4p; parameters <- newp; SE <- result_mcmc_newp$SD; ylim <- c(0,0.4)
    # result=resultNest_newp; p3p <- list(ylim=c(0, 0.5)); SE=result_mcmc_newp$SD
    # plotR(result=resultNest_newp, ylim=c(0, 0.5), SE=result_mcmc_newp$SD, xlimSE=c(25, 30))
    
    if (is.null(result) & is.null(c(parameters, fixed.parameters))) {
      stop("Or parameters or result from searchR must be provided !")
    }
    
    p3p <- list(...)
    # p3p <- list(NULL)
    
    # 3/4/2016
    SSM <- getFromNamespace(x=".SSM", ns="embryogrowth")
    afficheCI <- FALSE
    
    
    if (!is.list(col)) col <- list(col)
    if (!is.list(lty)) lty <- list(lty)
    if (!is.list(ltyCI)) ltyCI <- list(ltyCI)
    if (!is.list(lwd)) lwd <- list(lwd)
    if (!is.list(lwdCI)) lwdCI <- list(lwdCI)
    if (!is.list(set.par)) set.par <- list(set.par)
    if (!is.list(SE)) SE <- list(SE)
    if (!is.list(xlimR)) xlimR <- list(xlimR)
    if (!is.list(xlimSE)) xlimSE <- list(xlimSE)
    if (!is.list(legend)) legend <- list(legend)
    if (!is.list(parameters)) parameters <- list(parameters)
    if (!is.list(fixed.parameters)) fixed.parameters <- list(fixed.parameters)
    
    if (is.null(result)) result <- NA

    if (class(result)!="list") result <- list(result)
      
      nbr <- max(length(col), length(lty), length(ltyCI), length(lwd), length(lwdCI), 
                 length(set.par), length(SE), length(xlimR), length(xlimSE), 
                 length(legend), 
                 length(parameters), length(fixed.parameters), length(result)
                 )

      # Je les mets tous a la meme taille
      
      col <- as.list(rep(unlist(col), nbr)[1:nbr])
      lty <- as.list(rep(unlist(lty), nbr)[1:nbr])
      ltyCI <- as.list(rep(unlist(ltyCI), nbr)[1:nbr])
      lwd <- as.list(rep(unlist(lwd), nbr)[1:nbr])
      lwdCI <- as.list(rep(unlist(lwdCI), nbr)[1:nbr])
      set.par <- as.list(rep(unlist(set.par), nbr)[1:nbr])
      # je dois les grouper par deux
      xlimR <- rep(xlimR, nbr)[1:nbr]
      xlimSE <- rep(xlimSE, nbr)[1:nbr]
      
      
      SE <- c(SE, rep(NA, nbr-length(SE)))
      result <- c(result, rep(NA, nbr-length(result)))
      parameters <- c(parameters, rep(NA, nbr-length(parameters)))
      fixed.parameters <- c(fixed.parameters, rep(NA, nbr-length(fixed.parameters)))
      
      legend <- c(unlist(legend), rep("", nbr-length(unlist(legend))))
      
      for (i in 1:nbr) {
        if (length(result[[i]]) != 1) {
          if ((length(parameters[[i]]) <= 1) & ifelse(is.null(parameters[[i]]), TRUE, all(is.na(parameters[[i]])))) parameters[[i]] <- result[[i]]$par
          if ((length(fixed.parameters[[i]]) <= 1) & ifelse(is.null(fixed.parameters[[i]]), TRUE, all(is.na(fixed.parameters[[i]])))) fixed.parameters[[i]] <- result[[i]]$fixed.parameters
          if ((length(SE[[i]]) <= 1) & ifelse(is.null(SE[[i]]), TRUE, all(is.na(SE[[i]])))) SE[[i]] <- result[[i]]$SE
        }
      }

    
    premier <- TRUE
    output <- NULL
    
    # dans nbr j'ai le nombre de series
    for (rs in 1:nbr) {
      
      intermediaire <- data.frame(Temperatures=numeric(), Average=numeric(), CI.Minus=numeric(), CI.Plus=numeric())
      
      # J'introduis les paramtres fixes - 16/7/2012
        parssm <- c(parameters[[rs]], fixed.parameters[[rs]])
        model_p <- ifelse(all(names(parssm)!="Rho25") & all(names(parssm)!="k"), 
                          "anchor", "equation")

        # 4/4/2016 res peut etre une matrice si je suis en anchor
        res <- SE[[rs]]
        # je cree une variable pour savoir le type de modele
      
      # je suis en Anchor
      if (model_p == "anchor") {
        parvalue <- parssm[!(names(parssm) %in% c("rK", "K", "Scale"))]
        parvalueT <- as.numeric(names(parvalue))
        xlR <- xlimR[[rs]]
        if (any(is.na(xlR))) xlR <- c(min(parvalueT, na.rm=TRUE), max(parvalueT, na.rm=TRUE))
        if (xlR[1]>273) xlR <- xlR-273.15
        xlSE <- xlimSE[[rs]]
        if (any(is.na(xlSE))) xlSE <- c(min(parvalueT, na.rm=TRUE), max(parvalueT, na.rm=TRUE))
        if (xlSE[1]>273) xlSE <- xlSE-273.15
        
        
      } else {
        xlR <- xlimR[[rs]]
        xlSE <- xlimSE[[rs]]
      }
        
      
      if (length(xlR)==2) {
        x <- seq(from=xlR[1],to=xlR[2],length=100)
      } else {
        x <- xlR
      }
        
      # xSE <- seq(from=xlSE[1],to=xlSE[2],length=100)
      
      xSEDeb <- which.min(abs(x-xlSE[1]))
      xSEFin <- which.min(abs(x-xlSE[2]))
      xSE <- x[xSEDeb:xSEFin]
      
      
      # if (x<273) x <- x+273.15
      voutlist <- SSM(x, parssm)
      
      if (!is.null(size) & !is.na(parssm["transition_S"]) & !is.na(parssm["transition_P"])) {
        r <- voutlist[[1]]
        r_L <- voutlist[[2]]
        transition <- 1/(1+exp(parssm["transition_S"]*(size-parssm["transition_P"])))
        vout <- r*transition+r_L*(1-transition)
      } else {
        vout <- voutlist[[set.par[[rs]]]]
      }
      
      if (scaleY=="auto" & premier) scaleY <- 10^(-floor(log10(max(vout, na.rm = TRUE))))
      
      y <- vout
      
      
      
      if (premier) {
        L <- modifyList(list(type = "l", las=1, col=col[[rs]]
                             , lty=lty[[rs]], lwd=lwd[[rs]]
                             , axes = TRUE, bty = "n"
                             , xlab = expression("Temperatures in " * degree * "C")
                             , ylab = expression(paste0("r x", as.character(scaleY), "(mm."~min^-1~")")), xlim=xlim), 
                        modifyList(list(x=x, y=scaleY*y), p3p))
        do.call(plot, L) 
      } else {
        lines(x=x, y=scaleY*y, col=col[[rs]], lty=lty[[rs]], lwd=lwd[[rs]])
      }
      
      intermediaire <- rbind(intermediaire, data.frame(Temperature=x, Average=y, CI.Plus=rep(NA, length(y)), CI.Minus=rep(NA, length(y))))
            
###### J'affiche l'intervalle de confiance
      
      # J'affiche l'intervalle de confiance si je l'ai
      if (!is.null(res)) 
        if (!all(is.na(res))) {
        # Soit c'est directement une matrice qui est fournie
        if (class(res)=="matrix") {
        lines(xSE, SSM(xSE, res[1,])[[1]]*scaleY, type="l", col=col[[rs]], lty=ltyCI[[rs]], lwd=lwdCI[[rs]])
        lines(xSE, SSM(xSE, res[2,])[[2]]*scaleY, type="l", col=col[[rs]], lty=ltyCI[[rs]], lwd=lwdCI[[rs]])
        intermediaire[xSEDeb:xSEFin, "CI.Minus"] <- SSM(xSE, res[1,])[[1]]
        intermediaire[xSEDeb:xSEFin, "CI.Plus"] <- SSM(xSE, res[2,])[[2]]
        # Ou alors ce sont des parametres
      } else {
        # parssm c'est les parametres
        # res c'est les SD
        res2 <- parssm
        res2[] <- 0
        res2[match(names(res), names(res2))] <- res[na.omit(match(names(res2), names(res)))]
        res <- res2
        
        if (model_p=="anchor") {
          # Je suis en Anchor, je peux avoir des NA au début et à la fin
          
#          nm <- ! (names(res) %in% c("rK", "K", "Scale", "transition_S", "transition_P"))
#          newx <- res[nm]
#          newT <- as.numeric(names(res)[nm])
          
#          xSEp <- xSE[(xSE>min(newT[!is.na(newx)])-273.15) & (xSE<max(newT[!is.na(newx)])-273.15)]
          
          vout <- SSM(xSE, parssm[!is.na(res)]+x.SE*res[!is.na(res)])[[1]]
          
          lines(xSE, scaleY*vout, 
                type="l", col=col[[rs]], lty=ltyCI[[rs]], lwd=lwdCI[[rs]])
          
          intermediaire[xSEDeb:xSEFin, "CI.Plus"] <- vout
          
          vout <- SSM(xSE, parssm[!is.na(res)]-x.SE*res[!is.na(res)])[[1]]
          
          lines(xSE, scaleY*vout, 
                type="l", col=col[[rs]], lty=ltyCI[[rs]], lwd=lwdCI[[rs]])
          
          intermediaire[xSEDeb:xSEFin, "CI.Minus"] <- vout
          
        } else {
          
          if (!any(is.na(res))) {
          # Je suis en modèle SSM ou Weibull; il ne doit pas y avoir de NA

            # 8/2/2014 dans parssm j'ai les paramtres
            # Je cree une liste avec les parametres et la moyenne et moyenne^2 pour chaque x
            ess <- list(Parametre=matrix(rep(NA, length(parssm)*replicate.CI), ncol=length(parssm), dimnames=list(NULL, names(parssm))), moyenne=rep(0,length(xSE)), moyenne2=rep(0,length(xSE)))
            
            # Pour chacun des parametres, je tire la valeur dans une loi normale
            for (i in seq_along(parssm)) {
                ess$Parametre[,i] <- rnorm(replicate.CI, parssm[i], res[names(parssm[i])])
            }
            
            afficheCI <- TRUE
            
            # Pour chacun des replicats, je calcule l'attendu
            for (i in 1:replicate.CI) {
              valeurlist <- SSM(xSE+273.15, ess$Parametre[i,])
              
              if (!is.null(size) & !is.na(parssm["transition_S"]) & !is.na(parssm["transition_P"])) {
                r <- valeurlist[[1]]
                r_L <- valeurlist[[2]]
                transition <- 1/(1+exp(parssm["transition_S"]*(size-parssm["transition_P"])))
                valeur <- r*transition+r_L*(1-transition)
                
              } else {
                valeur <- valeurlist[[set.par[[rs]]]]
              }
              
              ess$moyenne <- ess$moyenne+valeur
              ess$moyenne2 <- ess$moyenne2+valeur^2
            }
            
            sdR=sqrt(ess$moyenne2/replicate.CI-(ess$moyenne/replicate.CI)^2)
            
            lines(xSE, (y-x.SE*sdR)*scaleY, type="l", col=col[[rs]], lty=ltyCI[[rs]], lwd=lwdCI[[rs]])
            lines(xSE, (y+x.SE*sdR)*scaleY, type="l", col=col[[rs]], lty=ltyCI[[rs]], lwd=lwdCI[[rs]])
            intermediaire[xSEDeb:xSEFin, "CI.Minus"] <- y-x.SE*sdR
            intermediaire[xSEDeb:xSEFin, "CI.Plus"] <- y+x.SE*sdR
# fin du test qu'il n'y a pas de NA dans le modele parametrique
          }
# fin du test SSM ou Weibull            
        }
# fin du test matrix de res
      }
# fin du test si res existe
      }

      par(new=TRUE)
      premier <- FALSE
      
      output <- c(output, list(intermediaire))
      
      # fin de la boucle des series de resultats
    }
    
    if (show.box) {
      if (afficheCI) {
        legend(local.box, c("Mean", "Confidence interval"), lty=c(lty[[1]], ltyCI[[1]]), lwd=c(lwd[[1]], lwdCI[[1]]), bty = "n")
      } else {
        legend(local.box, c("Mean"), lty=lty[[1]], lwd=lwd[[1]], bty = "n")
      }
    }
    
    
    if (show.legend & any(!is.na(unlist(legend)))) {
      legend("bottomright", unlist(legend), lty=unlist(lty), lwd=unlist(lwd), bty = "n", col=unlist(col))
    }
    
    if (show.anchors)
      if (all(names(parssm)!="Rho25") & all(names(parssm)!="k")) {
        parssm2 <- parssm[! (names(parssm) %in% c("rK", "K", "Scale", "transition_S", "transition_P"))]
        # je suis en anchor
        if (!is.null(pch.anchors) & !is.null(cex.anchors) & !is.null(col.anchors))
          points(x=as.numeric(names(parssm2))-273.15, y=rep(y.anchors, length(parssm2)), pch=pch.anchors, col=col.anchors, cex=cex.anchors) 
      }
    
    return(invisible(output))
  }
