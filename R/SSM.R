#########################################################
## r
## based on Schoolfield, R. M., Sharpe, P. J. & Magnuson, C. E. 
## 1981. Non-linear regression of biological temperature-dependent 
## rate models based on absolute reaction-rate theory. Journal of 
## Theoretical Biology, 88, 719-731.
##########################################################

.SSM <- function (T, parms, H=NULL) {
  
  # H est l'hygrométrie
  
  if (is.null(H)) H <- rep(0, length(T))
  
  # library("embryogrowth"); T <- c(25, 26, 30) ; parms <- resultNest_newp$par
  # getFromNamespace(".SSM", ns="embryogrowth")(c(25, 30), resultNest_newp$par)
  # plotR(resultNest_newp, ylim=c(0, 0.5))

  
  # Si je travaille en degree celsius, je convertis en Kelvin
  if (T[1]<273) T <- T + 273.15
  
  # print(parms)
  
  nm <- names(parms)
  
  if (any(nm=="Dallwitz_b1")) {
    
    if (is.na(parms["Dallwitz_b4"])) parms <- c(parms, Dallwitz_b4 = 6)
    if (is.na(parms["Dallwitz_b5"])) parms <- c(parms, Dallwitz_b5 =0.4)
    T <- T - 273.15
    c1 <- 1/(1+0.28*parms["Dallwitz_b4"]+0.72*log(1+parms["Dallwitz_b4"]))
    c2 <- 1+parms["Dallwitz_b4"]/(1+1.5*parms["Dallwitz_b4"]+0.39*parms["Dallwitz_b4"]^2)
    u <- (T-parms["Dallwitz_b3"])/(parms["Dallwitz_b3"]-parms["Dallwitz_b2"])-c1
    v <- (u + exp(parms["Dallwitz_b4"]*u))/c2
    rT <- 1E-5*abs(parms["Dallwitz_b1"]*10^(-v^2*(1-parms["Dallwitz_b5"]+parms["Dallwitz_b5"]*v^2)))
    rT_L <- rT
  } else {
  
  if (any(nm=="k")) {
    # C'est le modele de Weibull
    # Je suis en Kelvin
    rT <- dweibull(T, shape=abs(parms["k"]), scale=abs(parms["lambda"]))*parms["scale"]*1E-5
    if (any(nm=="k_L")) {
      rT_L <- dweibull(T, shape=abs(parms["k_L"]), scale=abs(parms["lambda_L"]))*parms["scale_L"]*1E-5
    } else {
      rT_L <- rT
    }
    
  } else 
    
    if (any(nm == "Peak")) {
      
      # Comme ca Peak est en degree Celsius
      T <- T - 273.15
      
      if (any(nm == "LengthB")) {
        
        xpar <- parms
        
        if (is.na(xpar["Flat"])) xpar["Flat"] <- 0
        
        if (!is.na(xpar["Length"])) xpar["LengthB"] <- xpar["LengthE"] <- xpar["Length"]
        
        xpar["MinB"] <- 0
        xpar["MinE"] <- 0
        xpar["Begin"] <- xpar["Peak"]-xpar["LengthB"]
        xpar["End"] <- xpar["Peak"]+xpar["LengthE"]
        
        xpar["PmoinsF"]<-xpar["Peak"]-(xpar["Flat"]/2)
        xpar["PplusF"]<-xpar["Peak"]+(xpar["Flat"]/2)
        
        xpar["PmoinsFB"]<-xpar["PmoinsF"]-xpar["Begin"]
        xpar["EPplusF"]<-xpar["End"]-xpar["PplusF"]
        
        xpar["MaxMinB"]<-xpar["Max"]-xpar["MinB"]
        xpar["MaxMinE"]<-xpar["Max"]-xpar["MinE"]
        
        # Modele sinusoidal
        rT <- 1E-5*ifelse(T<xpar["Begin"], xpar["MinB"],
                          ifelse(T<xpar["PmoinsF"], ((1+cos(pi*(xpar["PmoinsF"]-T)/xpar["PmoinsFB"]))/2)*xpar["MaxMinB"]+xpar["MinB"],
                                 ifelse(T<xpar["PplusF"], xpar["Max"],
                                        ifelse(T<xpar["End"], ((1+cos(pi*(T-(xpar["PplusF"]))/xpar["EPplusF"]))/2)*xpar["MaxMinE"]+xpar["MinE"],
                                               xpar["MinE"]
                                        )
                                 )
                          )
        )
      
} else {
      Tlogique <- ((T-parms["Peak"]) < 0)
      
      if (!is.na(parms["sd"])) parms["sdL"] <- parms["sdH"] <- parms["sd"]
      
      rT <- 1E-5*parms["Scale"]*(dnorm(T-parms["Peak"], mean=0, sd=ifelse(Tlogique, parms["sdL"], parms["sdH"]))*ifelse(Tlogique, 1, parms["sdH"]/parms["sdL"]))/dnorm(0, mean=0, sd=parms["sdL"])
    }
      
      rT_L <- rT
      
    } else {
    
    if (all(nm != "Rho25")) {
      # Je suis en Anchored
      # Je ne garde que les parametres necessaires
      newx <- parms[! (nm %in% c("rK", "K", "Scale", "transition_S", "transition_P"))]
      newx <- ifelse(newx<0, 0, newx)
      newT <- as.numeric(names(newx))
      scale <- ifelse(is.na(parms["Scale"]), 1, parms["Scale"])
      
      if (newT[1] < 273) newT <- 273.15+newT
      
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
      
     # listpolynom <- rep(list(NA), (length(newx)-3))
     # rangeX <- newT[2:(length(newx)-1)]
     # rangeX[1] <- rangeX[1]-10
     # rangeX[length(rangeX)] <- rangeX[length(rangeX)]+10
     # 
     #  for (i in seq_along(listpolynom)) {
     #    if (all(!is.na(newx[i:(i+3)]))) {
     #      listpolynom[[i]] <- list(
     #        polynom = poly.calc(newT[i:(i+3)], newx[i:(i+3)])
     #      )
     #    }
     #  }
     # 
     #  r <- rep(NA, length(T))
     #  pl <- findInterval(T, rangeX)
     # 
     #  for (i in seq_along(T)) {
     #    p <- listpolynom[[pl[i]]]$polynom
     #    if (!is.na(p[1])) r[i] <- predict(p, T[i])
     #  }
      
      
      rT <- ifelse(r<0, 1E-4, r)*scale*1E-5
      rT_L <- rT
      
    } else {
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
        rT<-( (rho25 + rho25_b*H) * (T/298)*exp((dha/R)*(unsur298-unsurT)))/(1+exp((dhh/R)*((1/t12H)-unsurT)))
      } else {
        # 28/7/2012 - T12H change en DT	
        #	"T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
        dhl <- parms["DHL"]*1E3
        t12L <- parms["T12L"]
        t12H <- t12L+abs(parms["DT"])
        rT<-((rho25 + rho25_b*H)*(T/298)*exp((dha/R)*(unsur298-unsurT)))/(1+exp((dhl/R)*((1/t12L)-unsurT))+exp((dhh/R)*((1/t12H)-unsurT)))
      }
      
      if (!is.na(parms["Rho25_L"])) {
        
        if (is.na(parms["Rho25_b_L"])) {rho25_b_L <- 0} else {rho25_b_L <- parms["Rho25_b_L"]/1E7}
        
        rho25_L <- parms["Rho25_L"]/1E7
        dha_L <- parms["DHA_L"]*1E3
        dhh_L <- parms["DHH_L"]*1E3
        
        if (is.na(parms["DHL_L"])) {
          #	"T12H", "DHA",  "DHH", "Rho25"
          t12H_L <- abs(parms["T12H_L"])
          rT_L <- ((rho25_L + rho25_b_L*H)*(T/298)*exp((dha_L/R)*(unsur298-unsurT)))/(1+exp((dhh_L/R)*((1/t12H_L)-unsurT)))
        } else {
          # 28/7/2012 - T12H change en DT	
          #	"T12L", "T12H", "DHA",  "DHH", "DHL", "Rho25"
          dhl_L <- parms["DHL_L"]*1E3
          t12L_L <- parms["T12L_L"]
          t12H_L <- t12L+abs(parms["DT_L"])
          rT_L <- ((rho25_L + rho25_b_L*H)*(T/298)*exp((dha_L/R)*(unsur298-unsurT)))/(1+exp((dhl_L/R)*((1/t12L_L)-unsurT))+exp((dhh_L/R)*((1/t12H_L)-unsurT)))
        }
      } else {
        rT_L <- rT
      }
    }
    }
  }
    
  epsilon <- parms["epsilon"]/1E7
  if (is.na(epsilon)) epsilon <- 0
  epsilon_L <- parms["epsilon_L"]/1E7
  if (is.na(epsilon_L)) epsilon_L <- 0
  
  return(list(unname(rT+epsilon), unname(rT_L+epsilon_L)))
  
}

.wrapperSSM <- function(...) {
  par <- unlist(list(...))
  return(getFromNamespace(".SSM", ns = "embryogrowth")(T=par["T"], parms=par[names(par) != "T"])[[1]])
}
