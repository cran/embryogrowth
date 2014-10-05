

.hessianRichardson<-function(x, pt) {
	return(numDeriv::hessian(.fonctionfit, x, method="Richardson", pt=pt))
	}
