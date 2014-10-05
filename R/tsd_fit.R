.tsd_fit <- function(par, males, N, temperatures, equation) {

if (equation=="logistic") {
  p <- 1/(1+exp((1/par["S"])*(par["P"]-temperatures)))
}

if (equation=="Richards") {
	p <- ifelse(par["K"]>3 & sign(par["P"]-temperatures)==sign(par["S"]), 
		0.5*exp((temperatures-par["P"])/(par["S"]*exp(par["K"]))), 
		(1+(2^exp(par["K"])-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(par["K"])))
}

if (equation=="Double-Richards") {
  p <- ifelse(temperatures<par["P"],
    ifelse(par["K1"]>3 & sign(par["P"]-temperatures)==sign(par["S"]), 
              0.5*exp((temperatures-par["P"])/(par["S"]*exp(par["K1"]))), 
              (1+(2^exp(par["K1"])-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(par["K1"]))),
    ifelse(par["K2"]>3 & sign(par["P"]-temperatures)==sign(par["S"]), 
           0.5*exp((temperatures-par["P"])/(par["S"]*exp(par["K2"]))), 
           (1+(2^exp(par["K2"])-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(par["K2"])))
  )
}


if (equation=="Hulin") {
  p <- (1+(2^exp(par["K1"]*temperatures+par["K2"])-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(par["K1"]*temperatures+par["K2"]))
}

if (equation=="Hill") {
  if (par["P"]>=0) {
  p <- 1/(1+exp((1/par["S"])*(log(par["P"])-log(temperatures))))
  } else {
    p <- NA
  }
}

  p <- ifelse(p==0, 1E-9, p)
  p <- ifelse(p==1, 1-1E-9, p)
    if (is.na(p[1])) {return(Inf)} else {
   return(-sum(dbinom(males, N, p, log = TRUE)))
    }
  
}
