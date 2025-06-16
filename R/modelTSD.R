.modelTSD <- function(par, temperatures, equation="logistic") {
  #  parx <- as.numeric(par)
  #  names(parx) <- colnames(par)
  #  par <- parx
  # embryogrowth:::.modelTSD(par, temperatures, equation)
  if (equation=="logistic")	{
    if (any(names(par)=="P_low")) {
      p <- 1/(1+exp((1/par["S_low"])*(par["P_low"]-temperatures))) * 1/(1+exp((1/par["S_high"])*(par["P_high"]-temperatures)))
    } else {
      p <- 1/(1+exp((1/par["S"])*(par["P"]-temperatures)))
    }
  }
  
  if (equation=="a-logistic") {
    
    if (2^exp(par["K"]) == 1) {
      p <- 1-0.5*(exp((1/par["S"])*(par["P"]-temperatures))^(-1/(exp(par["K"])-1)))
    } else {
      p <- (1+(2^exp(par["K"])-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(par["K"]))
    }
    
  }
  
  if (equation=="double-a-logistic") {
    lk1 <- log(2^exp(par["K1"])-1)
    if (is.infinite(lk1)) lk1 <- exp(par["K1"])*log(2)
    lk2 <- log(2^exp(par["K2"])-1)
    if (is.infinite(lk2)) lk2 <- exp(par["K2"])*log(2)
    
    llog1 <- log(1+exp(lk1 + (1/par["S"])*(par["P"]-temperatures)))
    llog1 <- ifelse(is.infinite(llog1), 
                    lk1 + (1/par["S"])*(par["P"]-temperatures), 
                    llog1)
    llog2 <- log(1+exp(lk2 + (1/par["S"])*(par["P"]-temperatures)))
    llog2 <- ifelse(is.infinite(llog2), 
                    lk2 + (1/par["S"])*(par["P"]-temperatures), 
                    llog2)
    
    p <- ifelse(temperatures<par["P"],
                exp((-1/exp(par["K1"]))*llog1) ,
                exp((-1/exp(par["K2"]))*llog2)
    )
  }
  
  if (equation=="logit")	{
    if (any(names(par)=="P_low")) {
      p <- 1/(1+exp(4*par["S_low"]*(par["P_low"]-temperatures))) * 1/(1+exp(4*par["S_high"]*(par["P_high"]-temperatures)))
    } else {
      p <- 1/(1+exp(4*par["S"]*(par["P"]-temperatures)))
    }
  }
  
  
  if (equation=="flexit") {
    if (any(names(par)=="P_low")) {
      par_low <- unname(par[c("P_low", "S_low", "K1_low", "K2_low")])
      names(par_low) <- c("P", "S", "K1", "K2")
      par_high <- unname(par[c("P_high", "S_high", "K1_high", "K2_high")])
      names(par_high) <- c("P", "S", "K1", "K2")
      p <- flexit(x=temperatures, par=par_low)*flexit(x=temperatures, par=par_high)
    } else {
      p <- flexit(x=temperatures, par=par)
    }
  }
  
  if (equation=="flexit*") {
    if (any(names(par)=="P_low")) {
      P_low <- par["P_low"]
      SL_low <- par["SL_low"]
      SH_low <- par["SH_low"]
      TransitionS_low <- par["TransitionS_low"]
      if (is.na(TransitionS_low)) TransitionS_low <- 100
      
      P_high <- par["P_high"]
      SL_high <- par["SL_high"]
      SH_high <- par["SH_high"]
      TransitionS_high <- par["TransitionS_high"]
      if (is.na(TransitionS_high)) TransitionS_high <- 100
      
      p <- 1/(1+exp(4*(SL_low /(1+exp(TransitionS_low*(temperatures - P_low))) + SH_low /(1+exp(TransitionS_low*(P_low - temperatures))))*(P_low-temperatures)))
      p <- p * 1/(1+exp(4*(SL_high /(1+exp(TransitionS_high*(temperatures - P_high))) + SH_high /(1+exp(TransitionS_high*(P_high - temperatures))))*(P_high-temperatures)))
    } else {
      P <- par["P"]
      SL <- par["SL"]
      SH <- par["SH"]
      TransitionS <- par["TransitionS"]
      if (is.na(TransitionS)) TransitionS <- 100
      p <- 1/(1+exp(4*(SL / (1+exp(TransitionS*(temperatures - P))) + SH / (1+exp(TransitionS*(P - temperatures))))*(P-temperatures)))
    }
  }
  
  if (equation=="stairs") {
    p <- ifelse(temperatures < par["TRTL"], 1, ifelse(temperatures > par["TRTH"], 0, invlogit(par["SRTRT"])))
  }
  
  if (equation=="flexit**") {
    if (any(names(par)=="P_low")) {
      P_low <- par["P_low"]
      if (!is.na(par["S_low"])) {
        SL_low <- par["S_low"]
        SH_low <- par["S_low"]
      } else {
        SL_low <- par["SL_low"]
        SH_low <- par["SH_low"]
      }
      
      TransitionS_low <- par["TransitionS_low"]
      l <- par["l"]
      if (is.na(TransitionS_low)) TransitionS_low <- 100
      if (is.na(l)) l <- 0.05
      QBT_low <- 1/ (1+ exp(TransitionS_low*(temperatures - P_low)))
      
      
      P_high <- par["P_high"]
      if (!is.na(par["S_high"])) {
        SL_high <- par["S_high"]
        SH_high <- par["S_high"]
      } else {
        SL_high <- par["SL_high"]
        SH_high <- par["SH_high"]
      }
      TransitionS_high <- par["TransitionS_high"]
      if (is.na(TransitionS_high)) TransitionS_high <- 100
      QBT_high <- 1/ (1+ exp(TransitionS_high*(temperatures - P_high)))
      
      p <- 1/(
        1+exp(
          ((-log((1-l)/l))/(SL_low*QBT_low+SH_low*(1-QBT_low)))*((P_low-temperatures))
        )
      )
      p <- p * 1/(
        1+exp(
          ((-log((1-l)/l))/(SL_high*QBT_high+SH_high*(1-QBT_high)))*((P_high-temperatures))
        )
      )
    } else {
      P <- par["P"]
      S <- par["S"]
      if (!is.na(par["S"])) {
        SL <- par["S"]
        SH <- par["S"]
      } else {
        SL <- par["SL"]
        SH <- par["SH"]
      }
      TransitionS <- par["TransitionS"]
      l <- par["l"]
      if (is.na(l)) l <- 0.05
      if (is.na(TransitionS)) TransitionS <- 100
      QBT <- 1/ (1+ exp(TransitionS*(temperatures - P)))
      
      p <- 1/(
        1+exp(
          ((-log((1-l)/l))/(SL*QBT+SH*(1-QBT)))*((P-temperatures))
        )
      )
    }
  }
  
  if (equation == "hill") if (par["P"] <= 0) {p <- rep(Inf, length(temperatures))} else {p <- 1/(1+exp((1/par["S"])*(log(par["P"])-log(temperatures))))}
  if (equation=="hulin") {
    Kx <- par["K1"]*temperatures+par["K2"]
    
    lk <- log(2^exp(Kx)-1)
    lk <- ifelse(is.infinite(lk), 
                 exp(Kx)*log(2), 
                 lk)
    llog <- log(1+exp(lk + (1/par["S"])*(par["P"]-temperatures)))
    llog <- ifelse(is.infinite(llog), 
                   lk + (1/par["S"])*(par["P"]-temperatures), 
                   llog)
    
    p <- exp((-1/exp(Kx))*llog)
  }
  if (equation=="gsd") p <- rep(0.5, length(temperatures))
  return(p)
  
}
