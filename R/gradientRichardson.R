
.gradientRichardson<-function(x, pt) {
	return(numDeriv::grad(.fonctionfit, x, method="Richardson", pt=pt))
	}
