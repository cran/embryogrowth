#########################################################
## r
## based on Schoolfield, R. M., Sharpe, P. J. & Magnuson, C. E. 
## 1981. Non-linear regression of biological temperature-dependent 
## rate models based on absolute reaction-rate theory. Journal of 
## Theoretical Biology, 88, 719-731.
##########################################################

.SSM <- function (T, parms, Hygrometry=NULL, verbose=FALSE) {
  
  # H est l'hygrométrie
  if (is.null(Hygrometry)) Hygrometry <- rep(0, length(T))
  
  # library("embryogrowth"); T <- c(25, 26, 30) ; parms <- resultNest_newp$par
  # getFromNamespace(".SSM", ns="embryogrowth")(c(25, 30), resultNest_newp$par)
  # plotR(resultNest_newp, ylim=c(0, 0.5))
  
  
  # Si je travaille en degree celsius, je convertis en Kelvin
  if (T[1] < 100) T <- T + 273.15
  
  rT_L <- rT <- rep(NA, length(T))
  
  # print(parms)
  
  nm <- names(parms)
  
  if (verbose) {
    print(d(parms))
  }
  
  if (any(nm=="Dallwitz_b1")) {
    
    if (is.na(parms["Dallwitz_b4"])) parms <- c(parms, Dallwitz_b4 = 6)
    if (is.na(parms["Dallwitz_b5"])) parms <- c(parms, Dallwitz_b5 =0.4)
    T_ec <- T - 273.15
    if (1+parms["Dallwitz_b4"] <= 0) {
      rT_L <- rT <- rep(NA, length(T_ec))
    } else {
      c1 <- 1/(1+0.28*parms["Dallwitz_b4"]+0.72*log(1+parms["Dallwitz_b4"]))
      
      
      c2 <- 1+parms["Dallwitz_b4"]/(1+1.5*parms["Dallwitz_b4"]+0.39*parms["Dallwitz_b4"]^2)
      u <- (T_ec-parms["Dallwitz_b3"])/(parms["Dallwitz_b3"]-parms["Dallwitz_b2"])-c1
      v <- (u + exp(parms["Dallwitz_b4"]*u))/c2
      rT <- 1E-5*abs(parms["Dallwitz_b1"]*10^(-v^2*(1-parms["Dallwitz_b5"]+parms["Dallwitz_b5"]*v^2)))
      rT_L <- rT
    }
    
  }
  
  if (any(nm=="k")) {
    # C'est le modele de Weibull
    # Je suis en Kelvin
    rT <- dweibull(T, shape=abs(parms["k"]), scale=abs(parms["lambda"]))*parms["scale"]*1E-5
    if (any(nm=="k_L")) {
      rT_L <- dweibull(T, shape=abs(parms["k_L"]), scale=abs(parms["lambda_L"]))*parms["scale_L"]*1E-5
    } else {
      rT_L <- rT
    }
    
  }
  
  if (any(nm == "Peak") & (any(nm == "LengthB") | any(nm == "Length"))) {
    
    # Comme ca Peak est en degree Celsius
    T_ec <- T - 273.15
    if (is.na(parms["Flat"])) parms["Flat"] <- 0
    
    if (!is.na(parms["Length"])) parms["LengthB"] <- parms["LengthE"] <- abs(parms["Length"])
    
    if (!is.na(parms["Min"]))      parms["MinB"] <- parms["MinE"] <- abs(parms["Min"])
    
    parms["LengthB"] <- abs(parms["LengthB"])
    parms["LengthE"] <- abs(parms["LengthE"])
    parms["MinB"] <- abs(parms["MinB"])
    parms["MinE"] <- abs(parms["MinE"])
    
    parms["Begin"] <- parms["Peak"]-parms["LengthB"]
    parms["End"] <- parms["Peak"]+parms["LengthE"]
    
    parms["PmoinsF"] <- parms["Peak"]-(parms["Flat"]/2)
    parms["PplusF"] <- parms["Peak"]+(parms["Flat"]/2)
    
    parms["PmoinsFB"] <- parms["PmoinsF"]-parms["Begin"]
    parms["EPplusF"] <- parms["End"]-parms["PplusF"]
    
    parms["MaxMinB"] <- parms["Max"]-parms["MinB"]
    parms["MaxMinE"] <- parms["Max"]-parms["MinE"]
    
    # Modele sinusoidal
    rT <- 1E-5*ifelse(T_ec<parms["Begin"], parms["MinB"],
                      ifelse(T_ec<parms["PmoinsF"], ((1+cos(pi*(parms["PmoinsF"]-T_ec)/parms["PmoinsFB"]))/2)*parms["MaxMinB"]+parms["MinB"],
                             ifelse(T_ec<parms["PplusF"], parms["Max"],
                                    ifelse(T_ec<parms["End"], ((1+cos(pi*(T_ec-(parms["PplusF"]))/parms["EPplusF"]))/2)*parms["MaxMinE"]+parms["MinE"],
                                           parms["MinE"]
                                    )
                             )
                      )
    )
    rT_L <- rT
  }
  
  if (any(nm == "Peak") & all(nm != "LengthB")) {
    Tlogique <- ((T_ec-parms["Peak"]) < 0)
    
    if (!is.na(parms["sd"])) parms["sdL"] <- parms["sdH"] <- parms["sd"]
    
    rT <- 1E-5*parms["Scale"]*(dnorm(T-parms["Peak"], mean=0, sd=ifelse(Tlogique, parms["sdL"], parms["sdH"]))*ifelse(Tlogique, 1, parms["sdH"]/parms["sdL"]))/dnorm(0, mean=0, sd=parms["sdL"])
    rT_L <- rT
  }
  
  
  if (length(na.omit(suppressWarnings(as.numeric(nm)))) != 0) {
    # Je suis en Anchored
    # Je ne garde que les parametres necessaires
    ppp <- suppressWarnings(as.numeric(nm))
    newx <- parms[!is.na(ppp)]
    newx <- ifelse(newx<0, 0, newx)
    newT <- as.numeric(names(newx))
    scale <- ifelse(is.na(parms["Scale"]), 1, parms["Scale"])
    
    if (newT[1] < 100) newT <- 273.15+newT
    
    # df <- data.frame(x=newT, y=newx)
    # # polyn <- lm(y ~ poly(x, degree = 2), data=df)
    # polyn <- loess(y ~x, data=df)
    # r <- predict(polyn, data.frame(x=T))
    
    newx <- c(newx[1], newx, tail(newx, n=1))
    newT <- c(newT[1]-10, newT, tail(newT, n=1)+10)
    
    listpolynom <- rep(list(NA), (length(newx)-1))
    rangeX <- newT
    
    for (i in seq_along(listpolynom)) {
      y <- newx[i:(i+1)]
      x <- newT[i:(i+1)]
      listpolynom[[i]] <- list(
        polynom = lm(y ~ x)
      )
    }
    
    r <- rep(NA, length(T))
    pl <- findInterval(T, rangeX)
    
    for (i in seq_along(T)) {
      p <- listpolynom[[pl[i]]]$polynom
      x <- T[i]
      r[i] <- predict(p, newdata=data.frame(x=T[i]))
    }
    
    
    rT <- ifelse(r<0, 1E-4, r)*scale*1E-5
    rT_L <- rT
    
  }
  
  if (any(nm == "Rho25")) {
    
    # je suis en model SSM
    
    if (is.na(parms["Rho25_b"])) {rho25_b <- 0} else {rho25_b <- parms["Rho25_b"]/1E7}
    
    R <- 8.314472
    rho25 <- parms["Rho25"]/1E7
    dha <- parms["DHA"]*1E3
    dhh <- parms["DHH"]*1E3
    unsurT <- 1/T
    unsur298 <- 0.00335570469798658 # 1/298
    
    if (is.na(parms["DHL"])) {
      # Je suis en 4 parametres
      #	"T12H", "DHA",  "DHH", "Rho25"
      t12H=abs(parms["T12H"])
      rT <- ( (rho25 + rho25_b*Hygrometry) * (T/298)*exp((dha/R)*(unsur298-unsurT)))/(1+exp((dhh/R)*((1/t12H)-unsurT)))
    } else {
      # 28/7/2012 - T12H change en DT	
      #	"T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
      dhl <- parms["DHL"]*1E3
      t12L <- parms["T12L"]
      t12H <- t12L+abs(parms["DT"])
      rT <- ((rho25 + rho25_b*Hygrometry)*(T/298)*exp((dha/R)*(unsur298-unsurT)))/(1+exp((dhl/R)*((1/t12L)-unsurT))+exp((dhh/R)*((1/t12H)-unsurT)))
    }
    
    if (!is.na(parms["Rho25_L"])) {
      
      if (is.na(parms["Rho25_b_L"])) {rho25_b_L <- 0} else {rho25_b_L <- parms["Rho25_b_L"]/1E7}
      
      rho25_L <- parms["Rho25_L"]/1E7
      dha_L <- parms["DHA_L"]*1E3
      dhh_L <- parms["DHH_L"]*1E3
      
      if (is.na(parms["DHL_L"])) {
        #	"T12H", "DHA",  "DHH", "Rho25"
        t12H_L <- abs(parms["T12H_L"])
        rT_L <- ((rho25_L + rho25_b_L*Hygrometry)*(T/298)*exp((dha_L/R)*(unsur298-unsurT)))/(1+exp((dhh_L/R)*((1/t12H_L)-unsurT)))
      } else {
        # 28/7/2012 - T12H change en DT	
        #	"T12L", "T12H", "DHA",  "DHH", "DHL", "Rho25"
        dhl_L <- parms["DHL_L"]*1E3
        t12L_L <- parms["T12L_L"]
        t12H_L <- t12L+abs(parms["DT_L"])
        rT_L <- ((rho25_L + rho25_b_L*Hygrometry)*(T/298)*exp((dha_L/R)*(unsur298-unsurT)))/(1+exp((dhl_L/R)*((1/t12L_L)-unsurT))+exp((dhh_L/R)*((1/t12H_L)-unsurT)))
      }
    } else {
      rT_L <- rT
    }
  }
  
  
  epsilon <- parms["epsilon"]/1E7
  if (is.na(epsilon)) epsilon <- 0
  epsilon_L <- parms["epsilon_L"]/1E7
  if (is.na(epsilon_L)) epsilon_L <- 0
  
  if (any(nm=="Threshold_Low")) {
    pth <- parms["Threshold_Low"]
    if (pth < 100) pth <- pth + 273.15
    rT <- ifelse(T < pth, 0, rT)
    rT_L <- ifelse(T < pth, 0, rT_L)
  }
  if (any(nm=="Threshold_High")) {
    pth <- parms["Threshold_High"]
    if (pth < 100) pth <- pth + 273.15
    rT <- ifelse(T > pth, 0, rT)
    rT_L <- ifelse(T > pth, 0, rT_L)
  }
  
  if (any(is.na(rT))) {
    stop(paste(c("Error in SSM function:\nT=", d(T), "\nparms=", d(parms)), collapse = " "))
  }
  
  return(list(unname(rT+epsilon), unname(rT_L+epsilon_L)))
  
}

# Attention hygrometry n'est pas pris en compte
.wrapperSSM <- function(...) {
  par <- list(...)
  set.par <- par[["set.par"]]
  if (is.null(set.par)) {
    set.par <- 1
  }
  Temp <- par[["T"]]
  if (is.null(Temp)) return(NA)
  parms <- par[(names(par) != "T") & (names(par) != "set.par")][[1]]
  return(getFromNamespace(".SSM", ns = "embryogrowth")(T=Temp, parms=parms)[[set.par]])
}
