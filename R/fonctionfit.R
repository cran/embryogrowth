.fonctionfit <- function(x, parameters=NULL, temperatures=NULL, derivate=NULL, 
                         weight=NULL,
                         test=NULL, M0=NULL, fixed.parameters=NULL) {

  
  return(info.nests(NestsResult=x, parameters=parameters, temperatures=temperatures,
                    derivate=derivate, weight=weight, test=test, 
                    M0=M0, fixed.parameters=fixed.parameters, 
                    out="Likelihood"))
}

