
.gradientRichardson <- function(x, temperatures, integral, weight, derivate, 
                              hatchling.metric, M0, fixed.parameters, WAIC=FALSE) {
	return(numDeriv::grad(info.nests, x, method="Richardson", temperatures=temperatures, 
	                      integral=integral, weight=weight, derivate=derivate, 
	                      hatchling.metric=hatchling.metric, M0=M0, fixed.parameters=fixed.parameters))
	}
