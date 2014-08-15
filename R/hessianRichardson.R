

.hessianRichardson<-function(x, pt) {
	return(hessian(.fonctionfit, x, method="Richardson", pt=pt))
	}
