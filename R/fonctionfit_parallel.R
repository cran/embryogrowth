.fonctionfit_parallel <- function(gl) {


# gl is list(y0=.EGR.env$nids[[namets]][1, "Mass"], K=Kval, R=r, R_L=r_L, 
#		transition_P=x[["transition_P"]], transition_S=x[["transition_S"]], series=namets))

Vfinale <- NULL

for (j in 1:length(gl)) {
#	namets <- gl[[j]]$series
	r <- gl[[j]]$R
	r_L <- gl[[j]]$R_L
	transition_P <- gl[[j]]$transition_P
	transition_S <- gl[[j]]$transition_S
	Kval <- gl[[j]]$K
	y <- gl[[j]]$y0
	tpk <- gl[[j]]$tempK
	tmin <- gl[[j]]$timemin
	der <- gl[[j]]$fnc
	meanSCL <- gl[[j]]$mean
	sdSCL <- gl[[j]]$sd
	

	for (i in 1:(length(tmin)-1)) {
			timesunique <- c(tmin[i], tmin[i+1])
			if (is.na(transition_S) | is.na(transition_P)) {
				transition <- 1
			} else {
				transition <- 1/(1+exp(transition_S*(y-transition_P)))
			}

			a <- as.numeric(r[tpk[i]]*transition+r_L[tpk[i]]*(1-transition))
			param <- c(alpha=a, K=Kval)
			out1 <- lsoda(y, timesunique, der, param)
			y <- as.numeric(out1[2,2])
		}
		
	L <- -dnorm(y, mean=meanSCL, sd=sdSCL, log=TRUE)
	names(L) <- gl[[j]]$series
	
	Vfinale <- c(Vfinale, L)
}

return(Vfinale)

}
