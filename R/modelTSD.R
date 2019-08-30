.modelTSD <- function(par, temperatures, equation="logistic") {
  #  parx <- as.numeric(par)
  #  names(parx) <- colnames(par)
  #  par <- parx
  # embryogrowth:::.modelTSD(par, temperatures, equation)
  if (equation=="logistic")	p <- 1/(1+exp((1/par["S"])*(par["P"]-temperatures)))
  
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
  
  if (equation=="flexit") {
    p <- flexit(x=temperatures, par=par)
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
