

.hessianRichardson<-function(x, temperatures, integral, derivate, weight,
                             hatchling.metric, M0, fixed.parameters) {
	return(numDeriv::hessian(info.nests, x, method="Richardson",
	                         temperatures=temperatures, 
	                         integral=integral, weight=weight, derivate=derivate, 
	                         hatchling.metric=hatchling.metric, M0=M0, fixed.parameters=fixed.parameters))
	}
