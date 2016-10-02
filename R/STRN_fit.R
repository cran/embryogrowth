.STRN_fit <- function(par, EmbryoGrowthTRN, tsd, Sexed, Males, Temperatures) {

  infoall <- info.nests(NestsResult=EmbryoGrowthTRN, SexualisationTRN=par, 
                        out="summary", replicate.CI = 1, progress=FALSE)$summary

    temp_TSD <- infoall[, Temperatures]
    names(temp_TSD) <- rev(rev(names(EmbryoGrowthTRN$data))[-(1:2)])
    
    # je dois les mettre dans le bon sens avec les noms
    Sexed <- Sexed[names(temp_TSD)]
    Males <- Males[names(temp_TSD)]
    
    
#    par_TSD <<- par
    
  -sum(dbinom(prob=predict(tsd, temp_TSD[!is.na(Sexed)], replicates=1)$sexratio, 
              size=Sexed[!is.na(Sexed)], x=Males[!is.na(Sexed)], log=TRUE))
    
}
