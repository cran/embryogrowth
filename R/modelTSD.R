.modelTSD <- function(par, temperatures, equation="logistic") {
#  parx <- as.numeric(par)
#  names(parx) <- colnames(par)
#  par <- parx
  # embryogrowth:::.modelTSD(par, temperatures, equation)
  if (equation=="logistic")	p <- 1/(1+exp((1/par["S"])*(par["P"]-temperatures)))
  if (equation=="richards") p <- ifelse(par["K"]>3 & sign(par["P"]-temperatures)==sign(par["S"]), 
                                        0.5*exp((temperatures-par["P"])/(par["S"]*exp(par["K"]))), 
                                        (1+(2^exp(par["K"])-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(par["K"])))
  if (equation=="double-richards") p <- ifelse(temperatures<par["P"], 
                                               ifelse(par["K1"]>3 & sign(par["P"]-temperatures)==sign(par["S"]), 
                                                      0.5*exp((temperatures-par["P"])/(par["S"]*exp(par["K1"]))), 
                                                      (1+(2^exp(par["K1"])-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(par["K1"]))) ,
                                               ifelse(par["K2"]>3 & sign(par["P"]-temperatures)==sign(par["S"]), 
                                                      0.5*exp((temperatures-par["P"])/(par["S"]*exp(par["K2"]))), 
                                                      (1+(2^exp(par["K2"])-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(par["K2"])))
  )
  if (equation=="hill") if (par["P"]<=0) {p <- rep(Inf, length(temperatures))} else {p <- 1/(1+exp((1/par["S"])*(log(par["P"])-log(temperatures))))}
  if (equation=="hulin") {
    Kx <- par["K1"]*temperatures+par["K2"]
    p <- ifelse(Kx>3 & sign(par["P"]-temperatures)==sign(par["S"]), 
                0.5*exp((temperatures-par["P"])/(par["S"]*exp(Kx))), 
                (1+(2^exp(Kx)-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(Kx)))
  }
  if (equation=="gsd") p <- rep(0.5, length(temperatures))
  return(p)
  
}
