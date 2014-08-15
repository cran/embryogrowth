#' tsd estimates the parameters that best describe temperature-dependent sex determination
#' @title Estimate the parameters that best describe temperature-dependent sex determination
#' @author Marc Girondot
#' @return A list the pivotal temperature, transitional range of temperatures and their SE
#' @param males A vector with male numbers
#' @param females A vector with female numbers
#' @param N A vector with total numbers
#' @param temperatures The constant incubation temperatures or any covariate used to fit sex ratio
#' @param df A dataframe with at least two columns named males, females or N and temperatures column
#' @param l The limit to define TRT (see Girondot, 1999)
#' @param parameters.initial Initial values for P, S or K search as a vector, ex. c(P=29, S=-0.3)
#' @param las.x las parameter for x axis
#' @param las.y las parameter for y axis
#' @param lab.PT Label to describe pivotal temperature
#' @param lab.TRT Label to describe transitional range of temperature
#' @param males.freq Should the graph uses males frequency [TRUE] or females [FALSE]
#' @param equation Could be "logistic", "Hill", "Richards" or "GSD"
#' @param replicate Number of replicate to estimate SE of TRT
#' @param range.CI The range of confidence interval for estimation, default=0.95
#' @param col.TRT The color of TRT
#' @param col.TRT.CI The color of CI of TRT based on range.CI
#' @param col.PT.CI The color of CI of PT based on range.CI
#' @param print Do the results must be printed at screen? TRUE or FALSE
#' @param ... Graphical parameters for plot(), exemple xlab="", ylab="", main=""
#' @description Estimate the parameters that best describe temperature-dependent sex determination
#' @references Girondot, M. 1999. Statistical description of temperature-dependent sex determination using maximum likelihood. Evolutionary Ecology Research, 1, 479-486.
#' @references Godfrey, M.H., Delmas, V., Girondot, M., 2003. Assessment of patterns of temperature-dependent sex determination using maximum likelihood model selection. Ecoscience 10, 265-272.
#' @references Hulin, V., Delmas, V., Girondot, M., Godfrey, M.H., Guillon, J.-M., 2009. Temperature-dependent sex determination and global change: are some species at greater risk? Oecologia 160, 493-506.
#' @examples
#' \dontrun{
#' CC_AtlanticSW <- subset(STSRE_TSD, RMU=="Atlantic, SW" & 
#'                           Species=="Caretta caretta" & Sexed!=0)
#' par(mar=c(4,4,5,1)+0.1)
#' tsdL <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="logistic"))
#' tsdH <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Hill"))
#' tsdR <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="Richards"))
#' gsd <- with (CC_AtlanticSW, tsd(males=Males, females=Females, 
#'                                  temperatures=Incubation.temperature-Correction.factor, 
#'                                  equation="GSD"))
#' compare_AIC(Logistic_Model=tsdL, Hill_model=tsdH, Richards_model=tsdR, GSD_model=gsd)
#' }
#' @export


tsd <- function(males=NULL, females=NULL, N=NULL, temperatures=NULL, 
	df=NULL, l=0.05, parameters.initial=c(P=NA, S=-0.5, K=0), males.freq=TRUE, 
	las.x=1, las.y=1, lab.PT="Pivotal temperature", 
	lab.TRT=paste0("Transitional range of temperatures l=",l*100,"%"), 
	col.TRT="gray", col.TRT.CI=rgb(0.8, 0.8, 0.8, 0.5), col.PT.CI=rgb(0.8, 0.8, 0.8, 0.5),
	equation="logistic", replicate=1000, range.CI=0.95, print=TRUE, ...) {
  
  # males <- c(10, 14, 7, 4, 3, 0, 0) ; females <- c(0, 1, 2, 4, 15, 10, 13)
  # temperatures <- c(25, 26, 27, 28, 29, 30, 31)
  # males <- NULL; females <- NULL; N <- NULL; temperatures <- NULL
  
  # N <- NULL; males <- CC_AtlanticSW$Males;females <- CC_AtlanticSW$Females; temperatures<- CC_AtlanticSW$Incubation.temperature-CC_AtlanticSW$Correction.factor
  # equation <- "Richards"
  # df <- NULL; l <- 0.05
  # parameters.initial=c(P=NA, S=-0.5, K=0); males.freq <- TRUE
  # las.x <- 1; las.y <- 1; lab.PT <- "Pivotal temperature"
  # lab.TRT <- "Transitional range of temperatures l=5%"; equation <- "logistic"
  # range.CI=0.95; replicate <- 1000; print=TRUE
  
  range.CI.qnorm <- qnorm(1-((1-range.CI)/2))
  
  if (!is.null(df)) {
    if (class(df)!="data.frame") {
      print("df parameter must be a data.frame")
      return(invisible())
    }
    
    namesdf <- tolower(names(df))
    males <- NULL
    females <- NULL
    N <- NULL
    
    if (any(namesdf=="males")) males <- df[, which(namesdf=="males")]
    if (any(namesdf=="male")) males <- df[, which(namesdf=="male")]
    if (any(namesdf=="females")) females <- df[, which(namesdf=="females")]
    if (any(namesdf=="female")) females <- df[, which(namesdf=="female")]
    if (any(namesdf=="n")) N <- df[, which(namesdf=="n")]
    if (any(namesdf=="temperatures")) temperatures <- df[, which(namesdf=="temperatures")]
    if (any(namesdf=="temperature")) temperatures <- df[, which(namesdf=="temperature")]
  }

  # si je retire des lignes, c'est décalé
#  if (!is.null(males)) males <- na.omit(males)
#  if (!is.null(females)) females <- na.omit(females)
#  if (!is.null(N)) N <- na.omit(N)
  if (is.null(temperatures)) {
    print("Temperatures must be provided")
    return(invisible())
#  } else {
    #    temperatures <- na.omit(temperatures)
  }
  
  cpt1 <- ifelse(is.null(males), 0, 1)+ifelse(is.null(females), 0, 1)+ifelse(is.null(N), 0, 1)
  if (cpt1<2) {
    print("At least two informations from males, females or N must be provided")
    return(invisible())
  }
  
  if (is.null(males)) males <- N-females
  if (is.null(females)) females <- N-males
  if (is.null(N)) N <- males+females
  
  if (length(temperatures)!=length(N)) {
    print("A temperature must be provided for each experiment")
    return(invisible())    
  }
  
  if (equation!="GSD" & equation!="logistic" & equation!="Hill" & equation!="Richards") {
    print("Equations supported are GSD, logistic, Hill and Richards")
    return(invisible())    
  }
  
  
  if (equation=="GSD") {
		result <- list(par = NULL, SE=NULL, hessian=NULL, 
		TRT=NULL,
		SE_TRT=NULL,
		message=NULL,
		convergence=0,
		counts=NULL,
		value=-sum(dbinom(x=males, size=males+females, prob=0.5, log=TRUE)),
		AIC=-2*sum(dbinom(x=males, size=males+females, prob=0.5, log=TRUE)))
		limit.low.TRT <- min(temperatures)
		limit.high.TRT <- max(temperatures)	
	} else {
    par <- parameters.initial
	
    if (is.na(par["P"])) {
      par["P"] <- temperatures[which.min(abs((males/(males+females))-0.5))]
    }
       if (equation=="logistic" | equation=="Hill") par <- par[which(names(par)!="K")]
    repeat {
      result  <- optim(par, .tsd_fit, males=males, N=N, temperatures=temperatures, equation=equation, method="BFGS", hessian=TRUE, control = list(maxit=1000))
      # result  <- optim(par, embryogrowth:::.tsd_fit, males=males, N=N, temperatures=temperatures, equation=equation, method="BFGS", hessian=TRUE, control = list(maxit=1000))

      if (result$convergence==0) break
      par<-result$par
      if (print) print("Convergence is not acheived. Optimization continues !")
    }
  
      par <- result$par
  
      mathessian <- result$hessian
  
      inversemathessian <- try(solve(mathessian), silent=TRUE)
  
  if (inherits(inversemathessian, "try-error")) {
    res <- rep(NA, length(par))
  } else {
    res <- abs(diag(inversemathessian))
  }
  
  names(res) <- names(par)
  result$SE <- res
  result$AIC <- 2*result$value+2*length(par)
  }
  
  if (equation!="GSD") {
	  l.l.TRT.c <- NULL
	  l.h.TRT.c <- NULL
	  TRT.c <- NULL
	
	for (j in 1:replicate) {
    
		par.P <- rnorm(1, par["P"], ifelse(j==1, 0, res["P"]))
    repeat {
		  par.S <- rnorm(1, par["S"], ifelse(j==1, 0, res["S"]))
      if (sign(par.S)==sign(par["S"])) break
    }
		if (equation=="Richards" & res["K"]!=0) {
		repeat {
		  par.K <- rnorm(1, par["K"], ifelse(j==1, 0, res["K"]))
		  if (sign(par.K)==sign(par["K"])) break
		}
    } else {par.K <- par["K"]}
    
		limit.low.TRT <- 20
		limit.high.TRT <- 40
		for (i in 1:5) {
	  		temperatures.se <- seq(from=limit.low.TRT, to=limit.high.TRT, length=i*20)
	  		if (equation=="logistic")	p <- 1/(1+exp((1/par.S)*(par.P-temperatures.se)))
 	  		if (equation=="Richards") p <- ifelse(par.K>3 & sign(par.P-temperatures.se)==sign(par.S), 
	  		                                      0.5*exp((temperatures.se-par.P)/(par.S*exp(par.K))), 
	  		                                      (1+(2^exp(par.K)-1)*exp((1/par.S)*(par.P-temperatures.se)))^(-1/exp(par.K)))
	  		if (equation=="Hill") p <- 1/(1+exp((1/par.S)*(log(par.P)-log(temperatures.se))))
			limit.low.TRT <- temperatures.se[tail(which(p>(1-l)), n=1)]
			limit.high.TRT <- temperatures.se[which(p < l)[1]]
			limit.low.TRT <- ifelse(identical(limit.low.TRT, numeric(0)), NA, limit.low.TRT)
			limit.high.TRT <- ifelse(identical(limit.high.TRT, numeric(0)), NA, limit.high.TRT)
      if (is.na(limit.low.TRT) | is.na(limit.high.TRT)) break
 		}
 		l.l.TRT.c <- c(l.l.TRT.c, limit.low.TRT)
 		l.h.TRT.c <- c(l.h.TRT.c, limit.high.TRT)
 		TRT.c <- c(TRT.c, limit.high.TRT-limit.low.TRT)
	}
  
	result$TRT <- as.numeric(TRT.c[1])
	result$TRT_limits <- as.numeric(c(l.l.TRT.c[1], l.h.TRT.c[1]))
	result$SE_TRT <- as.numeric(sd(TRT.c, na.rm = TRUE))
	result$SE_TRT_limits <- as.numeric(c(sd(l.l.TRT.c, na.rm = TRUE), sd(l.h.TRT.c, na.rm = TRUE)))
	limit.low.TRT <- result$TRT_limits[1]-range.CI.qnorm*result$SE_TRT_limits[1]
	limit.high.TRT <- result$TRT_limits[2]+range.CI.qnorm*result$SE_TRT_limits[2]
  }

   L <- list(...)

	if (males.freq) {
		L1 <- modifyList(list(x=temperatures, y=males/N, bty="n", type="n", xlab=expression("Temperatures in " * degree * "C"), ylab="male frequency"), L)
	} else {
		L1 <- modifyList(list(x=temperatures, y=females/N, bty="n", type="n", xlab=expression("Temperatures in " * degree * "C"), ylab="female frequency"), L)	
	}
  L1 <- modifyList(L1, list(ylim=c(0,1), xaxt="n", las=las.y))
  
  if (is.null(L$xlim)) {
    L1 <- modifyList(L1, list(xlim=c(floor(min(temperatures, limit.low.TRT)), floor(1+max(temperatures, limit.high.TRT)))))
  }
  
  a <- do.call(plot, L1) 
  
  x2 <- (par("usr")[1]+par("usr")[2]*26)/27
  x1 <- x2*26-par("usr")[2]/0.04
  
  axis(1, at=x1:x2, las=las.x)
  
  # je trace la TRT centrée sur P
  
  
if (equation!="GSD") {
  # limites de la TRT
  polygon(c(result$TRT_limits[1], result$TRT_limits[1], result$TRT_limits[2], result$TRT_limits[2]), c(0,1,1,0), border=NA, col=col.TRT)  
  # limites de la limite basse de la TRT
  polygon(c(result$TRT_limits[1]+range.CI.qnorm*result$SE_TRT_limits[1], result$TRT_limits[1]+range.CI.qnorm*result$SE_TRT_limits[1], result$TRT_limits[1]-range.CI.qnorm*result$SE_TRT_limits[1], result$TRT_limits[1]-range.CI.qnorm*result$SE_TRT_limits[1]), c(0,1,1,0), border=NA, col=col.TRT.CI)
  # limites de la limite haute de la TRT
  polygon(c(result$TRT_limits[2]+range.CI.qnorm*result$SE_TRT_limits[2], result$TRT_limits[2]+range.CI.qnorm*result$SE_TRT_limits[2], result$TRT_limits[2]-range.CI.qnorm*result$SE_TRT_limits[2], result$TRT_limits[2]-range.CI.qnorm*result$SE_TRT_limits[2]), c(0,1,1,0), border=NA, col=col.TRT.CI)
  # limites de la PT
  polygon(c(par["P"]-range.CI.qnorm*res["P"], par["P"]-range.CI.qnorm*res["P"], par["P"]+range.CI.qnorm*res["P"], par["P"]+range.CI.qnorm*res["P"]), c(0,1,1,0), border=NA, col=col.PT.CI)
  par(xpd=TRUE)
  segments(par["P"], 0, par["P"], 1.05, lty=4)
  segments(result$TRT_limits[1], 0, result$TRT_limits[1], 1.15, lty=3)
  segments(result$TRT_limits[2], 0, result$TRT_limits[2], 1.15, lty=3)
  text(x=par["P"], y=1.1, lab.PT)
  text(x=par["P"], y=1.2, lab.TRT)
  
}

  
	if (males.freq) {  
  		L1 <- modifyList(list(x=temperatures, y=males/N, bty="n", type="p", ylim=c(0,1)), L)
  	} else {
  		L1 <- modifyList(list(x=temperatures, y=females/N, bty="n", type="p", ylim=c(0,1)), L)
  	}
  L1 <- modifyList(L1, list(ylim=c(0,1), xlab="", ylab="", main="", axes=FALSE, xlim=c(x1, x2)))
  
  par(xpd=FALSE)

  par(new=TRUE)
  
  a <- do.call(plot, L1) 

	if (males.freq) {  
 	 	b <- binconf(males,N)
 	} else {
  		b <- binconf(females,N) 
	}
  # add error bars
  errbar(temperatures, b[,1], b[,3], b[,2], add=TRUE)
  
  
  mn <- min(temperatures)
  mx <- max(temperatures)
  
  par <- result$par
  
  x <- seq(from=x1, to=x2, length.out=100)
	if (equation=="logistic") {
		p <-  1/(1+exp((1/par["S"])*(par["P"]-x)))
		prop.theorique <- 1/(1+exp((1/par["S"])*(par["P"]-temperatures)))
	}
	if (equation=="Hill") {
		p <-  1/(1+exp((1/par["S"])*(log(par["P"])-log(x))))
		prop.theorique <- 1/(1+exp((1/par["S"])*(log(par["P"])-log(temperatures))))
	}
	if (equation=="GSD") {
		p <-  rep(0.5,length(x))
		prop.theorique <- rep(0.5,length(temperatures))
	}
	if (equation=="Richards") {
	p <- ifelse(par["K"]>3 & sign(par["P"]-x)==sign(par["S"]), 
			0.5*exp((x-par["P"])/(par["S"]*exp(par["K"]))), 
			(1+(2^exp(par["K"])-1)*exp((1/par["S"])*(par["P"]-x)))^(-1/exp(par["K"])))
	prop.theorique <- ifelse(par["K"]>3 & sign(par["P"]-temperatures)==sign(par["S"]), 
			0.5*exp((temperatures-par["P"])/(par["S"]*exp(par["K"]))), 
			(1+(2^exp(par["K"])-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(par["K"])))
	}

 # prop.cumul <- NULL
 # for (i in 1:length(N)) prop.cumul <- c(prop.cumul, binom.test(males[i], N[i], prop.theorique[i])$p.value)
 # 
 # 
 # result$GOF <- combine.test(prop.cumul, method = c("z.transform"))

	saturated <- 0
  for (i in 1:length(N)) saturated <- saturated- dbinom(males[i], N[i], males[i]/N[i], log=TRUE)
  LRT <- 2*(result$value-saturated)
  result$GOF <- 1-pchisq(LRT, length(N)-length(par))
  
  par(new=TRUE)
	if (males.freq) {   
  L1 <- modifyList(list(x=x, y=p, bty="n"), L)
  } else {
  L1 <- modifyList(list(x=x, y=1-p, bty="n"), L)
  }
  L1 <- modifyList(L1, list(ylim=c(0,1), axes=FALSE, xlab="", ylab="", type="l", main="", xlim=c(x1, x2)))

  a <- do.call(plot, L1)
  
  	pm <- NULL
		pp <- NULL

  
	if (equation=="logistic") {
    # ici je laisse S contrairement à predict.tsd
	  parameter.expand <- expand.grid(S=c(ifelse(sign(par["S"]-range.CI.qnorm*res["S"])!=sign(par["S"]), sign(par["S"])*1e-5, par["S"]-range.CI.qnorm*res["S"]), 
	                                      ifelse(sign(par["S"]+range.CI.qnorm*res["S"])!=sign(par["S"]), sign(par["S"])*1e-5, par["S"]+range.CI.qnorm*res["S"])),
	                                  P=c(par["P"]-range.CI.qnorm*res["P"], par["P"]+range.CI.qnorm*res["P"]), temperatures=x)
	  y <- as.numeric(with(parameter.expand, 1/(1+exp((1/S)*(P-temperatures)))))
	  l <- split(y, parameter.expand$temperatures)
	  
	  pm <- unlist(lapply(l, min))
	  pp <- unlist(lapply(l, max))
	}
  
	if (equation=="Hill") {
	  # ici je laisse S contrairement à predict.tsd
    parameter.expand <- expand.grid(S=c(ifelse(sign(par["S"]-range.CI.qnorm*res["S"])!=sign(par["S"]), sign(par["S"])*1e-5, par["S"]-range.CI.qnorm*res["S"]), 
                                        ifelse(sign(par["S"]+range.CI.qnorm*res["S"])!=sign(par["S"]), sign(par["S"])*1e-5, par["S"]+range.CI.qnorm*res["S"])),
                                    P=c(par["P"]-range.CI.qnorm*res["P"], par["P"]+range.CI.qnorm*res["P"]), temperatures=x)
    y <- as.numeric(with(parameter.expand, 1/(1+exp((1/S)*(log(P)-log(temperatures))))))
    l <- split(y, parameter.expand$temperatures)
    
    pm <- unlist(lapply(l, min))
    pp <- unlist(lapply(l, max))
	}
	
	if (equation=="Richards") {
	  parameter.expand <- expand.grid(S=c(ifelse(sign(par["S"]-range.CI.qnorm*res["S"])!=sign(par["S"]), sign(par["S"])*1e-5, par["S"]-range.CI.qnorm*res["S"]), 
	                                  ifelse(sign(par["S"]+range.CI.qnorm*res["S"])!=sign(par["S"]), sign(par["S"])*1e-5, par["S"]+range.CI.qnorm*res["S"])),
	                                  P=c(par["P"]-range.CI.qnorm*res["P"], par["P"]+range.CI.qnorm*res["P"]), 
                                    K=c(ifelse(sign(par["K"]-range.CI.qnorm*res["K"])!=sign(par["K"]), sign(par["K"])*1e-5, par["K"]-range.CI.qnorm*res["K"]), 
                                        ifelse(sign(par["K"]+range.CI.qnorm*res["K"])!=sign(par["K"]), sign(par["K"])*1e-5, par["K"]+range.CI.qnorm*res["K"])), 
	                                  temperatures=x)
	  y <- as.numeric(with(parameter.expand, 
	                       ifelse(K>3 & sign(P-temperatures)==sign(S), 
	                              0.5*exp((temperatures-P)/(S*exp(K))), 
	                              (1+(2^exp(K)-1)*exp((1/S)*(P-temperatures)))^(-1/exp(K)))
                         ))
	  l <- split(y, parameter.expand$temperatures)
	  
	  pm <- unlist(lapply(l, min))
	  pp <- unlist(lapply(l, max))

	}
	
  
  
	if (!is.null(pp)) {
	par(new=TRUE)
  	if (males.freq) {   

  L1 <- modifyList(list(x=x, y=pm, bty="n"), L)
  } else {
  L1 <- modifyList(list(x=x, y=1-pm, bty="n"), L)
  
  }
  L1 <- modifyList(L1, list(ylim=c(0,1), axes=FALSE, xlab="", ylab="", type="l", main="", lty=2, xlim=c(x1, x2)))
  a <- do.call(plot, L1) 
  
  par(new=TRUE)
  	if (males.freq) {   

  L1 <- modifyList(list(x=x, y=pp, bty="n"), L)
  } else {
  L1 <- modifyList(list(x=x, y=1-pp, bty="n"), L)
  
  }
  L1 <- modifyList(L1, list(ylim=c(0,1), axes=FALSE, xlab="", ylab="", type="l", main="", lty=2, xlim=c(x1, x2)))
  a <- do.call(plot, L1) 

	}  
  
  
  result$males <- males
  result$females <- females
  result$N <- N
  result$temperatures <- temperatures
  result$males.freq <- males.freq
  result$equation <- equation
  result$l <- l
  
  if (print) print(paste("The goodness of fit test is", sprintf("%.3f",result$GOF)))

if (equation!="GSD" & print) {
  print(paste("The", lab.PT, "is", sprintf("%.3f",par["P"]), "SE", sprintf("%.3f",res["P"])))
  print(paste("The", lab.TRT, "is", sprintf("%.3f",result$TRT), "SE", sprintf("%.3f",result$SE_TRT)))
  print(paste("The lower limit of ", lab.TRT, "is", sprintf("%.3f",result$TRT_limits[1]), "SE", sprintf("%.3f",result$SE_TRT_limits[1])))
  print(paste("The higher limit of ", lab.TRT, "is", sprintf("%.3f",result$TRT_limits[2]), "SE", sprintf("%.3f",result$SE_TRT_limits[2])))
}

	class(result) <- "tsd"
  return(invisible(result))
  
  
}
