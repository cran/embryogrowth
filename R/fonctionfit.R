.fonctionfit <- function(x, pt) {

# 12/11/2012 rajout du test pour un seul processeur fictif si pas parallel
if (pt$parallel) {
	nc <- detectCores()
} else {
	nc <- 1
}

# dans x j'ai les paramètres à ajuster
# Il faut que je rajoute les fixe - 16/7/2012
x <- c(x, pt$fixed.parameters)

NbTS <- pt$temperatures[["IndiceT"]][3]

tempK <- as.numeric(pt$temperatures[["Temperatures"]])
rlist <- .SSM(tempK, x)
r <- rlist[[1]]
r_L <- rlist[[2]]

names(r)<-tempK

names(r_L)<-tempK


gl <- NULL

p <- NULL
if (!is.na(x["K"])) {
	Kval <- x[["K"]]
	} else {
	Kval <- NULL
}




for (j in 1:NbTS) {
	namets <- names(pt$temperatures)[j]

	if (!is.na(x["rK"])) {
		Kval <- x[["rK"]]*pt$test[namets, "Mean"]
	}
	

	gl <- c(gl, list(list(y0=pt$temperatures[[namets]][1, "Mass"], K=Kval, R=r, R_L=r_L, 
		transition_P=x["transition_P"][[1]], transition_S=x["transition_S"][[1]],
		series=namets, tempK= as.character(pt$temperatures[[namets]][, "Temperatures K"]), 
			timemin=pt$temperatures[[namets]][,"Time"], fnc= pt$derivate, 
			mean=as.numeric(pt$test[namets, "Mean"]), sd=as.numeric(pt$test[namets, "SD"]))))

}

# nombre de série ns=lenght(gl)
# nombre de groupes=nc
# donc floor(ns/nc) serie par groupe et jusqu'au dernier pour le dernier

# ng, nombre de séries par groupe
ng <- floor(length(gl)/nc)

# 12/11/2012 rajout du test pour un seul processeur

if (ng!=0 & nc!=1) {

for (j in 1:(nc-1)) {

glencours <- gl[((j-1)*ng+1):(j*ng)]

if (pt$parallel) {
	s <- mcparallel(.fonctionfit_parallel(glencours))
} else {
	s <- .fonctionfit_parallel(glencours)
}

p <- c(p, list(s))

}
}

# je fais le dernier

j <- nc
glencours <- gl[((j-1)*ng+1):length(gl)]

if (pt$parallel) {
	s <- mcparallel(.fonctionfit_parallel(glencours))
} else {
	s <- .fonctionfit_parallel(glencours)
}

p <- c(p, list(s))
 
if (pt$parallel) {
	p <- mccollect(p)
	names(p) <- NULL
}

L <- unlist(p)
# dans L j'ai un vecteur avec le nom
# il faut que j'applique le pt$weight

if (!is.null(pt$weight)) {
	L <- L[order(names(L))]*pt$weight[order(names(pt$weight))]
}


return(sum(L))



}
