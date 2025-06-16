#' FormatNests creates a dataset of class "Nests" to be used with searchR
#' @title Create a dataset of class Nests to be used with searchR
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with all the nests formated to be used with searchR.
#' @param data Data to be newly formated.
#' @param Time.Format Format of time. See description. If NULL, no time conversion is done.
#' @param Time.Zone The format of time zone as obtained by OlsonNames(). See description.
#' @param previous Data already formated.
#' @param LayingTime Named POSIXct or POSIXlt time for each nest in data.
#' @param UnitTime The units for time as a named list or vector
#' @param Longitude The longitude of the nests as a named list or vector
#' @param Latitude The latitude of the nests as a named list or vector
#' @param Informations Some textual information about the nests as a named list or vector
#' @param Males Number of sexed eggs being males
#' @param Females Number of sexed eggs being females
#' @param simplify If TRUE, simply the time series by removing identical time series of temperatures.
#' @param usemiddletime If TRUE, suppose that recorded temperatures are those at middle segment.
#' @param weight Named vector with weight used to estimate likelihood.
#' @param hatchling.metric.mean The average size of hatchlings
#' @param hatchling.metric.sd The standard deviation of size of hatchlings
#' @param col.Time Name of the column with time.
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
#' The parameter Weight is a vector: weight=c(Nest1=1, Nest2=1.2).\cr
#' The parameter LayingTime is also a vector of POSIXct time or POSIXlt time.\cr
#' It can be used to format database already formated with old format; in this case, 
#' just use data=xxx with xxx being the old format database.\cr
#' The UnitTime should be "seconds", "minutes", "hours", or "days" to be understood by plot function.\cr
#' @examples 
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(data=nest, previous=NULL, col.Time="Time")
#' # If I try to add the same nest, I have an error
#' formated <- FormatNests(data=nest, previous=formated, col.Time="Time")
#' # I duplicate the database and change the names
#' nest_duplicate <- nest
#' colnames(nest_duplicate) <- paste0(colnames(nest_duplicate), "_essai")
#' formated <- FormatNests(data=nest_duplicate, previous=formated, col.Time="Time_essai")
#' # It is possible to add information about these nests
#' formated <- FormatNests(data=nest, previous=NULL, col.Time="Time")
#' formated <- UpdateNests(data=formated, Males=c(DY.1=10), Females=c(DY.1=2))
#' ####################
#' Laying.Time <- matrix(c("DY.1", "15/05/2010", 
#'                  "DY.17", "24/05/2010", 
#'                  "DY.16", "24/05/2010", 
#'                  "DY.18", "25/05/2010", 
#'                  "DY.20", "25/05/2010", 
#'                  "DY.21", "26/05/2010", 
#'                  "DY.22", "26/05/2010", 
#'                  "DY.23", "26/05/2010", 
#'                  "DY.24", "27/05/2010", 
#'                  "DY.25", "27/05/2010", 
#'                  "DY.28", "28/05/2010", 
#'                  "DY.26", "28/05/2010", 
#'                  "DY.27", "28/05/2010", 
#'                  "DY.146", "20/06/2010", 
#'                  "DY.147", "20/06/2010", 
#'                  "DY.172", "24/06/2010", 
#'                  "DY.175", "24/06/2010", 
#'                  "DY.170", "24/06/2010", 
#'                  "DY.260", "06/07/2010", 
#'                  "DY.282", "12/07/2010", 
#'                  "DY.310", "18/07/2010", 
#'                  "DY.309", "18/07/2010", 
#'                  "DY.328", "25/07/2010", 
#'                  "DY.331", "26/07/2010"), byrow=TRUE, ncol=2)
#' tz <- OlsonNames()[grepl("Asia/Istanbul", OlsonNames())]
#' Laying.Time_f <- setNames(as.POSIXlt.character(Laying.Time[, 2], format = "%d/%m/%Y", tz=tz), 
#'                            Laying.Time[, 1])
#' formated <- FormatNests(data=nest, previous=NULL, col.Time="Time", LayingTime=Laying.Time_f)
#' ####################
#' # Now when the data are with absolute dates that are already formatted
#' nest_ec <- data.frame(Time=as.POSIXlt("24/05/2010", format="%d/%m/%Y")+ nest[, 1]*60, 
#'                       DY.1.x=nest[, 2])
#' formated <- FormatNests(data=nest_ec, previous=NULL, col.Time="Time")
#' ####################
#' # Now when the data are with absolute date that are in text format for example after 
#' #   reading a csv format
#' nest_ec <- data.frame(Time=format(as.POSIXlt("24/05/2010", format="%d/%m/%Y")+ nest[, 1]*60, 
#'                                   format = "%d/%m/%Y %H:%M:%S"), 
#'                       DY.1.x=nest[, 2])
#' formated <- FormatNests(data=nest_ec, previous=NULL, col.Time="Time", 
#'                         Time.Format="%d/%m/%Y %H:%M:%S", 
#'                         Time.Zone=OlsonNames()[grepl("Asia/Istanbul", OlsonNames())], 
#'                         hatchling.metric.mean=39.33, hatchling.metric.sd=1.92)
#' 
#' }
#' @export


FormatNests <-
  function(data=stop("A dataset must be provided !"), 
           Time.Format = NULL                       ,
           Time.Zone=NULL                           ,
           previous=NULL                            , 
           LayingTime=NULL                          ,
           UnitTime="minutes"                       ,
           Longitude=NULL                           ,
           Latitude=NULL                            ,
           Informations=NULL                        ,
           Males=NULL                               ,
           Females=NULL                             ,
           usemiddletime=FALSE                      , 
           simplify=TRUE                            , 
           weight=NULL                              , 
           hatchling.metric.mean=NULL               ,
           hatchling.metric.sd=NULL                 ,
           col.Time="Time"                          ) {
    
    # Je cree une une fonction qui prepare un fichier pour etre utilise
    # Les differents nids seront des matrices dans une liste
    # Dans chaque matrice on a le temps depuis le debut de l'incubation
    # la temperature en C, la temperature en K, la valeur de r et la masse
    
    # previous=NULL; simplify=TRUE; weight=NULL
    
    
    UnitTime <- match.arg(tolower(UnitTime), choices = c("minutes", "hours", "seconds"))
    
    if (inherits(previous, "Nests")) { 
      previous <- UpdateNests(previous) 
    }
    
    if (inherits(data, "Nests")) { 
      data <- UpdateNests(data)
    }
    
    
    # Si data est déjà un format Nest2, je ne fais rien
    if (!inherits(data, "Nests2")) {
      
      nests_f <- list()
      # meme chose que nids=as.list(NULL)
      
      names_nests <- colnames(data)[colnames(data) != col.Time]
      Time <- data[, col.Time]
      if (!is.null(Time.Format)) {
        Time <- as.POSIXlt(data[, col.Time], format=Time.Format, tz=Time.Zone)
      }
      
      if (length(unique(names_nests)) != length(names_nests)) {
        stop("All nest names must be different; check the names.")
      }
      
      for (j in names_nests) {
        
        # je veux le transformer en une liste avec a l'interieur un dataframe par enregistreur
        # je commence par sortir la liste un par un qui n'ont pas de NA
        
        # je stocke les donnees de l'enregistreur dans une matrice
        
        newess2 <- as.numeric(subset(data[,j], !is.na(data[,j])))
        newessT <- Time[!is.na(data[,j])]
        
        LayingTime_ec <- LayingTime[j]
        
        if (is.null(LayingTime_ec)) {
          LayingTime_ec <- Time[1]
        }
        
        if (is.numeric(LayingTime_ec)) {
          LayingTime_ec <- NULL
        } else {
          if (!is.numeric(newessT[1])) {
            if (UnitTime=="seconds") newessT <- as.numeric(newessT) - as.numeric(newessT[1])
            if (UnitTime=="minutes") newessT <- as.numeric(newessT)/60 - as.numeric(newessT[1])/60
            if (UnitTime=="hours") newessT <- as.numeric(newessT)/(60*60) - as.numeric(newessT[1])/(60*60)
          }
        }
        
        # 19/10/2012 je calcule les etats intermediaires en terme de temps
        
        ess2 <- newess2[1]
        essT <- newessT[1]
        
        for(i in 1:(length(newess2)-1)) {
          ess2 <- c(ess2, newess2[i+1])
          if (usemiddletime) {
            minutes <- newessT[i]+(newessT[i+1]-newessT[i])/2
          } else {
            minutes <- newessT[i+1]
          }
          essT <- c(essT, minutes)
          
        }
        
        ess2 <- c(ess2, newess2[length(newess2)])
        essT <- c(essT, newessT[length(newess2)])
        
        ess <- matrix(c(essT, ess2, ess2+273.15, rep(NA, 3*length(ess2))), ncol=6)
        
        colnames(ess)<-c("Time", "Temperatures C", "Temperatures K", "r", "Mass", "IndiceK")
        
        
        # ensuite je supprime les temps avec des valeurs de temperatures identiques
        # sauf la derniere qui doit rester - 20/7/2012
        # si que deux lignes, je ne fais rien - 27/7/2012
        if ((dim(ess)[1]>2) & simplify) {
          for (i in 2:(dim(ess)[1]-1)) {
            if (ess[i,2]==ess[i-1,2]) ess[i,1]=NA
          }
          
          if (!is.na(ess[nrow(ess)-1,1]))
            if (ess[nrow(ess), 1] == ess[nrow(ess)-1,1]) {
              ess[nrow(ess),1]=NA
            }
        }
        
        # nidsEC[[names(data[j])]]<-subset(ess, !is.na(ess[,1]))
        # 23/4/2015
        nests_f_ec <- list()
        nests_f_ec$data <- subset(ess, !is.na(ess[,1]))
        if (!is.null(names(weight))) {
          nests_f_ec$weight <- weight[j]
        } else {
          nests_f_ec$weight <- 1
        }
        nests_f_ec$LayingTime <- LayingTime_ec
        nests_f_ec$UnitTime <- UnitTime
        if (!is.null(names(Longitude))) {
          if (!(j %in% names(Longitude))) 
            warning(paste0("The timeseries of name ", j, " is not among the values for Longitude."))
          
          nests_f_ec$Longitude <- Longitude[j]
        } else {
          nests_f_ec$Longitude <- Longitude[1]
        }
        if (!is.null(names(Latitude))) {
          if (!(j %in% names(Latitude))) 
            warning(paste0("The timeseries of name ", j, " is not among the values for Latitude."))
          
          nests_f_ec$Latitude <- Latitude[j]
        } else {
          nests_f_ec$Latitude <- Latitude[1]
        }
        if (!is.null(names(Informations))) {
          if (!(j %in% names(Informations))) 
            warning(paste0("The timeseries of name ", j, " is not among the values for Informations."))
          
          nests_f_ec$Informations <- Informations[j]
        } else {
          nests_f_ec$Informations <- Informations
        }
        if (!is.null(names(Males))) {
          if (!(j %in% names(Males))) 
            warning(paste0("The timeseries of name ", j, " is not among the values for Males."))
          
          nests_f_ec$Males <- Males[j]
        } else {
          nests_f_ec$Males <- Males
        }
        if (!is.null(names(Females))) {
          if (!(j %in% names(Females))) 
            warning(paste0("The timeseries of name ", j, " is not among the values for Females."))
          
          nests_f_ec$Females <- Females[j]
        } else {
          nests_f_ec$Females <- Females
        }
        if (!is.null(names(hatchling.metric.mean))) {
          if (!(j %in% names(hatchling.metric.mean))) 
            warning(paste0("The timeseries of name ", j, " is not among the values for hatchling.metric.mean."))
          nests_f_ec$hatchling.metric.mean <- hatchling.metric.mean[j]
        } else {
          nests_f_ec$hatchling.metric.mean <- hatchling.metric.mean[1]
        }
        if (!is.null(names(hatchling.metric.sd))) {
          if (!(j %in% names(hatchling.metric.sd))) 
            warning(paste0("The timeseries of name ", j, " is not among the values for hatchling.metric.sd."))
          nests_f_ec$hatchling.metric.sd <- hatchling.metric.sd[j]
        } else {
          nests_f_ec$hatchling.metric.sd <- hatchling.metric.sd[1]
        }
        nests_f_ec$Name <- j
        nests_f_ec$Temperatures <- unique(as.character(nests_f_ec$data[, "Temperatures K"]))
        nests_f_ec$IndiceT <- c(Tmin=min(nests_f_ec$data[, "Temperatures K"]), Tmax=max(nests_f_ec$data[, "Temperatures K"]), NbTS=1)
        nests_f_ec <- list(Nests = nests_f_ec)
        names(nests_f_ec) <- j
        nests_f <- c(nests_f, nests_f_ec)
      }
      
      data <- list(Nests=nests_f, 
                   Names=names(nests_f), 
                   Temperatures=unique(unlist(lapply(nests_f, FUN=function(x) return(x$Temperatures)))), 
                   IndiceT=c(Tmin=min(unlist(lapply(nests_f, FUN=function(x) return(x$IndiceT["Tmin"])))), 
                             Tmax=max(unlist(lapply(nests_f, FUN=function(x) return(x$IndiceT["Tmax"])))), 
                             NbTS=length(nests_f))
      )
    }
    
    if (!is.null(previous)) {
      
      if (length(unique(c(data$Names, previous$Names))) != length(c(data$Names, previous$Names))) {
        stop("All nest names must be different; check the names.")
      }
      
      nests_f <- c(data$Nests, previous$Nests)
      
      data <- list(Nests=nests_f, 
                   Names=names(nests_f), 
                   Temperatures=unique(unlist(lapply(nests_f, FUN=function(x) return(x$Temperatures)))), 
                   IndiceT=c(Tmin=min(unlist(lapply(nests_f, FUN=function(x) return(x$IndiceT["Tmin"])))), 
                             Tmax=max(unlist(lapply(nests_f, FUN=function(x) return(x$IndiceT["Tmax"])))), 
                             NbTS=length(nests_f))
      )
    }
    
    data <- addS3Class(data, "Nests2")
    
    # Je cherche toutes les températures
    # Je l'ai c'est dans data$Temperatures
    temp <- data$Temperatures
    
    for (j in data$Names) {
      data$Nests[[j]]$data[, "IndiceK"] <- match(as.character(data$Nests[[j]]$data[, "Temperatures K"]), temp)
    }
    
    return(data)
  }
