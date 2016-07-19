

#########################################################
## r
## based on Schoolfield, R. M., Sharpe, P. J. & Magnuson, C. E. 
## 1981. Non-linear regression of biological temperature-dependent 
## rate models based on absolute reaction-rate theory. Journal of 
## Theoretical Biology, 88, 719-731.
##########################################################

.SSM <- function (T, parms) {
  
  # library("embryogrowth"); T <- c(25, 26, 30) ; parms <- resultNest_newp$par
  # getFromNamespace(".SSM", ns="embryogrowth")(c(25, 30), resultNest_newp$par)
  # plotR(resultNest_newp, ylim=c(0, 0.5))

  
  # Si je travaille en degree celsius, je convertis en Kelvin
  if (T[1]<273) T <- T+273.15
  
  nm <- names(parms)
  
  if (any(nm=="k")) {
    rT <- dweibull(T-parms["theta"], shape=abs(parms["k"]), scale=abs(parms["lambda"]))*parms["scale"]*1E-5
    rT_L <- rT
    
  } else
    
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
        rT<-(rho25*(T/298)*exp((dha/R)*(unsur298-unsurT)))/(1+exp((dhh/R)*((1/t12H)-unsurT)))
      } else {
        # 28/7/2012 - T12H change en DT	
        #	"T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
        dhl <- parms["DHL"]*1E3
        t12L <- parms["T12L"]
        t12H <- t12L+abs(parms["DT"])
        rT<-(rho25*(T/298)*exp((dha/R)*(unsur298-unsurT)))/(1+exp((dhl/R)*((1/t12L)-unsurT))+exp((dhh/R)*((1/t12H)-unsurT)))
      }
      
      if (!is.na(parms["Rho25_L"])) {
        
        rho25_L <- parms["Rho25_L"]/1E7
        dha_L <- parms["DHA_L"]*1E3
        dhh_L <- parms["DHH_L"]*1E3
        
        if (is.na(parms["DHL_L"])) {
          #	"T12H", "DHA",  "DHH", "Rho25"
          t12H_L <- abs(parms["T12H_L"])
          rT_L <- (rho25_L*(T/298)*exp((dha_L/R)*(unsur298-unsurT)))/(1+exp((dhh_L/R)*((1/t12H_L)-unsurT)))
        } else {
          # 28/7/2012 - T12H change en DT	
          #	"T12L", "T12H", "DHA",  "DHH", "DHL", "Rho25"
          dhl_L <- parms["DHL_L"]*1E3
          t12L_L <- parms["T12L_L"]
          t12H_L <- t12L+abs(parms["DT_L"])
          rT_L <- (rho25_L*(T/298)*exp((dha_L/R)*(unsur298-unsurT)))/(1+exp((dhl_L/R)*((1/t12L_L)-unsurT))+exp((dhh_L/R)*((1/t12H_L)-unsurT)))
        }
      } else {
        rT_L <- rT
      }
    }
  
  return(list(as.numeric(rT), as.numeric(rT_L)))
  
}
