.STRN_fit <- function(par, fixed.parameters=NULL, 
                      equation=equation, 
                      TSP.borders, embryo.stages, 
                      EmbryoGrowthTRN, tsd, Sexed, Males, Temperatures, 
                      zero = 1E-9, 
                      parallel) {
  
  serafaire <- names(Sexed[(!is.na(Sexed)) & !is.na(Males)])
  
  parSTRN <- c(par, fixed.parameters)
  parSTRN <- parSTRN[!(names(parSTRN) %in% c("P", "S", "K", "K1", "K2"))]
  
  if (identical(parSTRN, structure(numeric(0), .Names = character(0)))) {
    parSTRN <- NULL
  }
  
  infoall <- info.nests(NestsResult=EmbryoGrowthTRN, 
                        series=serafaire, 
                        embryo.stages=embryo.stages, 
                        TSP.borders=TSP.borders, 
                        SexualisationTRN=parSTRN, 
                        out="summary", replicate.CI = 0, 
                        progress=FALSE, 
                        parallel=parallel)$summary
  
  temp_TSD <- infoall[, Temperatures]
  names(temp_TSD) <- rownames(infoall)
  
  # je dois les mettre dans le bon sens avec les noms
  Sexed <- Sexed[serafaire]
  Males <- Males[serafaire]
  temp_TSD <- temp_TSD[serafaire]
  
  if (is.null(tsd)) {
    
    tsd <- list(par = c(par, fixed.parameters), equation=equation, type="temperature")
    class( tsd) <- "tsd"
    
  }
  sr <- predict(tsd, temperatures = temp_TSD, replicate.CI = NULL, probs=0.5)[1, ]
  sr <- ifelse(sr==0, zero, sr)
  sr <- ifelse(sr==1, 1-zero, sr)

  return(-sum(dbinom(prob=sr, 
              size=Sexed, x=Males, log=TRUE))
  )
}
