#' UpdateNests creates a dataset of class "Nests2" to be used with searchR
#' @title Create a dataset of class Nests2 from an object of class Nests
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with all the nests formated to be used with searchR.
#' @param data An object of class Nests or Nests2.
#' @param LayingTime Named POSIXct or POSIXlt time for each nest in data.
#' @param UnitTime The units for time as a named list or vector
#' @param Longitude The longitude of the nests as a named list or vector
#' @param Latitude The latitude of the nests as a named list or vector
#' @param Informations Some textual information about the nests as a named list or vector
#' @param Males Number of sexed eggs being males.
#' @param Females Number of sexed eggs being females.
#' @param hatchling.metric.mean The average size of hatchlings.
#' @param hatchling.metric.sd The standard deviation of size of hatchlings.
#' @param weight The weight of different nests for likelihood estimation.
#' @description Will create a dataset of class Nests2 to be used with searchR\cr
#' This function is used to convert Nests or Nests2 format into the new one, Nests2, and add information.\cr
#' @examples 
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' nest <- FormatNests(nest)
#' nest2 <- UpdateNests(data=nest)
#' nest2 <- UpdateNests(data=nest2, Males=c(DY.1=10), Females=c(DY.1=20))
#' }
#' @export


UpdateNests <-
  function(data=stop("An object with nests must be provided !"), 
           weight=NULL                              ,
           LayingTime=NULL                          ,
           UnitTime=NULL                            ,
           Longitude=NULL                           ,
           Latitude=NULL                            ,
           Informations=NULL                        ,
           Males=NULL                               ,
           Females=NULL                             , 
           hatchling.metric.mean=NULL               ,
           hatchling.metric.sd=NULL                 ) {
    
    if (!inherits(data, "Nests")) { 
      
      if (!inherits(data, "Nests2")) { 
      stop("The object is not of required format.")
      } else {
        # Je suis déjà en Nests2
        nests_names <- data$Names
        
        for (i in nests_names) {
          data[["Nests"]][[i]][["weight"]] <- modifyVector(data[["Nests"]][[i]][["weight"]], weight[i])
          data[["Nests"]][[i]][["LayingTime"]] <- modifyVector(data[["Nests"]][[i]][["LayingTime"]], LayingTime[i])
          data[["Nests"]][[i]][["UnitTime"]] <- modifyVector(data[["Nests"]][[i]][["UnitTime"]], UnitTime[i])
          data[["Nests"]][[i]][["Longitude"]] <- modifyVector(data[["Nests"]][[i]][["Longitude"]], Longitude[i])
          data[["Nests"]][[i]][["Latitude"]] <- modifyVector(data[["Nests"]][[i]][["Latitude"]], Latitude[i])
          data[["Nests"]][[i]][["Informations"]] <- modifyVector(data[["Nests"]][[i]][["Informations"]], Informations[i])
          data[["Nests"]][[i]][["Males"]] <- modifyVector(data[["Nests"]][[i]][["Males"]], Males[i])
          data[["Nests"]][[i]][["Females"]] <- modifyVector(data[["Nests"]][[i]][["Females"]], Females[i])
        }
        return(data)
        
      }
      
    } else {
      weight <- modifyVector(data[["weight"]], weight)
      LayingTime <- modifyVector(data[["LayingTime"]], LayingTime)
      UnitTime <- modifyVector(data[["UnitTime"]], UnitTime)
      Longitude <- modifyVector(data[["Longitude"]], Longitude)
      Latitude <- modifyVector(data[["Latitude"]], Latitude)
      Informations <- modifyVector(data[["Informations"]], Informations)
      Males <- modifyVector(data[["Males"]], Males)
      Females <- modifyVector(data[["Females"]], Females)
      
      nests_names <- names(data)
      nests_names <- nests_names[!(nests_names %in% c("weight1", "LayingTime", "UnitTime", "Longitude", "Latitude", "Informations", "IndiceT", "Temperatures"))]
      nests_f <- list()
      for(i in nests_names) {
        nests_f_ec <- list()
        # names(nests_f_ec) <- i
        nests_f_ec$data <- data[[i]]
        nests_f_ec$weight <- data[["weight"]][i]
        nests_f_ec$LayingTime <- data[["LayingTime"]][i]
        nests_f_ec$UnitTime <- data[["UnitTime"]][i]
        nests_f_ec$Longitude <- data[["Longitude"]][i]
        nests_f_ec$Latitude <- data[["Latitude"]][i]
        nests_f_ec$Informations <- data[["Informations"]][i]
        nests_f_ec$Males <- data[["Males"]][i]
        nests_f_ec$Females <- data[["Females"]][i]
        nests_f_ec$hatchling.metric.mean <- hatchling.metric.mean[i]
        nests_f_ec$hatchling.metric.sd <- hatchling.metric.sd[i]
        nests_f_ec$Name <- i
        nests_f_ec$Temperatures <- unique(as.character(data[[i]][, "Temperatures K"]))
        nests_f_ec$IndiceT <- c(Tmin=min(data[[i]][, "Temperatures K"]), Tmax=max(data[[i]][, "Temperatures K"]), NbTS=1)
        nests_f_ec <- list(Nests = nests_f_ec)
        names(nests_f_ec) <- i
        nests_f <- c(nests_f, nests_f_ec)
      }
      
      nests_f <- list(Nests=nests_f, 
                      Names=names(nests_f), 
                      Temperatures=unique(unlist(lapply(nests_f, FUN=function(x) return(x$Temperatures)))), 
                      IndiceT=c(Tmin=min(unlist(lapply(nests_f, FUN=function(x) return(x$IndiceT["Tmin"])))), 
                                Tmax=max(unlist(lapply(nests_f, FUN=function(x) return(x$IndiceT["Tmax"])))), 
                                NbTS=length(nests_f))
      )
      
      nests_f <- addS3Class(nests_f, "Nests2")
    }
    return(nests_f)
    }
