# deeply modified from a MCMC script from Olivier Martin (INRA, Paris-Grignon)

# Algo Metropolis-Hastings
# ------------------------

.MHalgoGen<-function(n.iter=10000, parameters=NULL, data=NULL, likelihood=NULL, 
n.chains = 4, n.adapt = 100, thin=30, trace=FALSE)
{

  require("coda")
  
t <- as.character(trace)
pt <- NULL
if (t=="TRUE") {pt <- 1;tf <- TRUE}
if (t=="FALSE") {pt <- 0;tf <- FALSE}
if (is.null(pt)) {
  tf <- TRUE
  pt <- floor((n.adapt+n.iter)/trace)
}

cpt_trace <- 0


res<-as.list(NULL)
resL<-as.list(NULL)

for (kk in 1:n.chains) {

# Initialisation
nbvar<-dim(parameters)[1]
varp<-matrix(rep(NA, (nbvar+1)*(n.adapt+n.iter+2)), ncol=nbvar+1)
varp2<-matrix(rep(NA, (nbvar+1)*(n.adapt+n.iter+2)), ncol=nbvar+1)

colnames(varp)<-c(rownames(parameters), "Ln L")
colnames(varp2)<-c(rownames(parameters), "Ln L")


# for(i in 1:nbvar) varp[1, i]<-as.numeric(parameters[i, 'Init'])
varp[1, 1:nbvar]<-as.numeric(parameters[1:nbvar, 'Init'])

varp[1, "Ln L"]<- (-likelihood(data, varp[1, 1:nbvar]))
cpt<-1
varp2[cpt, 1:nbvar]<-varp[1, 1:nbvar]
varp2[cpt, "Ln L"]<-varp[1, "Ln L"]
cpt<-2


if (trace) {
	cat(paste("Chain ", kk, ": [", 1, "] ",as.numeric(varp[1, nbvar+1]), "\n", sep=""))
} else {
	cat(paste("Chain ", kk, "\n", sep=""))
}

MaxL<-varp[1,]

sdg=NULL

# 18/1/2013
# for(i in 1:nbvar) sdg<-c(sdg, as.numeric(parameters[i, 'SDProp']))
sdg<-c(sdg, as.numeric(parameters[1:nbvar, 'SDProp']))

Prior<-matrix(as.numeric(parameters[,2:3]), ncol=2)

Limites<-matrix(as.numeric(parameters[,5:6]), ncol=2)

dfun<-parameters[,"Density"]


# Itérations
for (i in 2:(n.adapt+n.iter+1))
{
	newvarp<-varp[i-1, 1:nbvar]
	for (j in 1:nbvar) {	
		# Nouveau paramètre
		propvarp<-newvarp
		propvarp[j]<-propvarp[j]+rnorm(1, mean=0, sd=sdg[j])

		if (propvarp[j]<=Limites[j,2] && propvarp[j]>=Limites[j,1]) 
			{
			logratio<-(get(dfun[j])(propvarp[j],Prior[j,1],Prior[j,2],log=T)+
					-likelihood(data, propvarp)-
		        	(get(dfun[j])(newvarp[j],Prior[j,1],Prior[j,2],log=T)+varp[i-1, "Ln L"]))
			alpha<-min(c(1,exp(logratio)))
			if (runif(1, min=0, max=1)<=alpha) {newvarp<-propvarp} 
			}
	}		
	varp[i, 1:nbvar]<-newvarp
	varp[i, "Ln L"]<-(-likelihood(data, newvarp))
	
	if (MaxL["Ln L"]<varp[i, "Ln L"]) {MaxL<-varp[i,]}
  

# 6/10/2012: Je stocke tout	
#	if ((i>n.adapt) && ((i%%thin)==0)) {
		varp2[cpt, 1:nbvar]<-newvarp
		varp2[cpt, "Ln L"]<-varp[i, "Ln L"]
		cpt<-cpt+1
#	}
	if (tf) {
    cpt_trace <- cpt_trace+1
    if (cpt_trace>=pt) {
	    cat(paste("Chain ", kk, ": [", i, "] ", as.numeric(varp[i, "Ln L"]), "\n", sep=""))
      cpt_trace <- 0
    }
  }
}

lp <- coda::as.mcmc(varp2[1:(cpt-1), 1:nbvar])
lp <- coda::mcmc(data=lp, start=n.adapt+1, end=n.iter, thin=thin)

res<-c(res, list(lp))
resL <- c(resL, list(varp2[1:(cpt-1), "Ln L"]))

}


names(res) <- 1:n.chains

res <- coda::as.mcmc.list(res)

cat("Best likelihood for: \n")
for (j in 1:nbvar) {cat(paste(names(MaxL[j]), "=", MaxL[j], "\n"))}



out <- (list(resultMCMC=res, resultLnL=resL, parametersMCMC=list(parameters=parameters, n.iter=n.iter, n.chains=n.chains, n.adapt=n.adapt, thin=thin)))
class(out) <- "mcmcComposite"
return(out)

}