#' FormatNests creates a dataset of class "Nests" to be used with searchR
#' @title Create a dataset of class Nests to be used with searchR
#' @author Marc Girondot
#' @return A list with all the nests formated to be used with searchR.
#' @param data Data to be newly formated
#' @param previous Data already formated
#' @param simplify If TRUE, simply the time series by removing identical time series of temperatures
#' @param weight Named vector with weight for likelihhod
#' @description Will create a dataset of class Nests to be used with searchR\cr
#' FormatNests(nest, previous=x) with x being a previously formated data.\cr
#' The raw data must be organized being:\cr
#' First column is the time in minutes since the beginning of incubation\cr
#' Each column next is the trace of temperatures, one column for each nest.\cr
#' For example, for two nests:\cr
#' Time   Nest1    Nest2\cr
#' 0       29.8     27.6\cr
#' 90      30.2     28.8\cr
#' 120     30.4     30.7\cr
#' 180     31.2     32.6\cr
#' ...\cr
#' 65800   30.8     32.6\cr
#' 65890            30.2\cr
#' 65950            30.4\cr
#' \cr
#' The Nest1 ends incubation at 65800 minutes whereas Nest2 ends incubation at 65950 (last row\cr
#' with temperature for each).\cr
#' The parameter Weight is a vector: weight=c(Nest1=1, Nest2=1.2)\cr
#' \cr
#' @examples 
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest, previous=NULL)
#' formated <- FormatNests(nest)
#' }
#' @export


FormatNests <-
function(data=stop("A dataset must be provided !"), previous=NULL, simplify=TRUE, weight=NULL) {

# Je crée une fonction qui prépare un fichier pour être utilisé
# Les différents nids seront des matrices dans une liste
# Dans chaque matrice on a le temps depuis le début de l'incubation
# la tepérature en °C, la température en K, la valeur de r et la masse

# previous=NULL; simplify=TRUE; weight=NULL

nidsEC <- vector(mode="list", length=0)
# même chose que nids=as.list(NULL)

for (j in 2:dim(data)[2]) {

# je veux le transformer en une liste avec à l'intérieur un dataframe par enregistreur
# je commence par sortir la liste un par un qui n'ont pas de NA

# je stocke les données de l'enregistreur dans une matrice

newess2 <- as.numeric(subset(data[,j], !is.na(data[,j])))
newessT <- as.numeric(subset(data[,1], !is.na(data[,j])))

# 19/10/2012 je calcule les états intermédaires en terme de temps

ess2 <- newess2[1]
essT <- newessT[1]

for(i in 1:(length(newess2)-1)) {
	ess2 <- c(ess2, newess2[i+1])
	minutes <- newessT[i]+(newessT[i+1]-newessT[i])/2
	essT <- c(essT, minutes)
	
}

ess2 <- c(ess2, newess2[length(newess2)])
essT <- c(essT, newessT[length(newess2)])

ess<-matrix(c(essT, ess2, ess2+273.15, rep(NA, 2*length(ess2))), ncol=5)

colnames(ess)<-c("Time", "Temperatures C", "Temperatures K", "r", "Mass")


# ensuite je supprime les temps avec des valeurs de températures identiques
# sauf la dernière qui doit rester - 20/7/2012
# si que deux lignes, je ne fais rien - 27/7/2012
if ((dim(ess)[1]>2) & simplify) {
	for (i in 2:(dim(ess)[1]-1)) {
    	if (ess[i,2]==ess[i-1,2]) ess[i,1]=NA
	}
}

nidsEC[[names(data[j])]]<-subset(ess, !is.na(ess[,1]))

}


if (!is.null(previous)) {
	nidsprevious<-previous[1:previous$IndiceT["NbTS"]]
	nidsEC<-c(nidsEC, nidsprevious)
}

# J'enregistre les tempmin et max car toujours pareil - Je n'en ai plus besoin; je les garde cependant

# j'enregestre les temp comme des facteurs
nbts <- length(nidsEC)
temp<-NULL
for (j in 1:nbts) {
	temp<- c(temp, nidsEC[[j]][, "Temperatures K"])
}
tempmin <- min(temp)
tempmax <- max(temp)

tempaf<- as.factor(temp)
templevels <- levels(tempaf)

nidsEC[["IndiceT"]] <- c(Tmin=tempmin,Tmax=tempmax, NbTS=nbts)
nidsEC[["Temperatures"]] <- templevels

nidsEC[["weight"]] <- weight

class(nidsEC) <- "Nests"

return(nidsEC)

}
