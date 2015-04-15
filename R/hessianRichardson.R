

.hessianRichardson<-function(x, temperatures, derivate, weight,
                             test, M0, fixed.parameters) {
	return(numDeriv::hessian(.fonctionfit, x, method="Richardson",
	                         temperatures=temperatures, 
	                         derivate=derivate, weight=weight,
	                         test=test, M0=M0, fixed.parameters=fixed.parameters))
	}
