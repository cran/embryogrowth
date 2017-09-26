.STRN_fit <- function(par, fixed.parameters=NULL, EmbryoGrowthTRN, tsd, Sexed, Males, Temperatures, parallel) {
  
  serafaire <- names(Sexed[!is.na(Sexed)])
  
  infoall <- info.nests(NestsResult=EmbryoGrowthTRN, series=serafaire, 
                        SexualisationTRN=c(par, fixed.parameters), 
                        out="summary", replicate.CI = 1, progress=FALSE, 
                        parallel=parallel)$summary
  
  infoall2 <- info.nests(NestsResult=EmbryoGrowthTRN, series=serafaire, 
                        SexualisationTRN=c(par, fixed.parameters), 
                        out="summary", replicate.CI = 1, progress=FALSE, 
                        parallel=parallel)$summary
  
  temp_TSD <- infoall[, Temperatures]
  names(temp_TSD) <- rownames(infoall)
  
  # je dois les mettre dans le bon sens avec les noms
  Sexed <- Sexed[serafaire]
  Males <- Males[serafaire]
  temp_TSD <- temp_TSD[serafaire]

  -sum(dbinom(prob=predict(tsd, temperatures=temp_TSD, replicate.CI=NULL, probs=0.5)[1, ], 
              size=Sexed, x=Males, log=TRUE))
  
}
