.tsd_fit <- function(par, males, N, temperatures, equation) {

if (equation=="logistic") {
  p <- 1/(1+exp((1/par["S"])*(par["P"]-temperatures)))
}

if (equation=="Richards") {
	p <- ifelse(par["K"]>3 & sign(par["P"]-temperatures)==sign(par["S"]), 
		0.5*exp((temperatures-par["P"])/(par["S"]*exp(par["K"]))), 
		(1+(2^exp(par["K"])-1)*exp((1/par["S"])*(par["P"]-temperatures)))^(-1/exp(par["K"])))
}

if (equation=="Hill") {
  p <- 1/(1+exp((1/par["S"])*(log(par["P"])-log(temperatures))))
}

  p <- ifelse(p==0, 1E-9, p)
  p <- ifelse(p==1, 1-1E-9, p)

   return(-sum(dbinom(males, N, p, log = TRUE)))
  
}
