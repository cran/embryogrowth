.fonctionSTRNMCMC <- function(data, x) {

# .STRN_fit <- function(par, EmbryoGrowthTRN, tsd, Sexed, Males, Temperatures) {
  
return(getFromNamespace(".STRN_fit", ns="embryogrowth")(par=x, 
                 fixed.parameters=data$fixed.parameters, 
                 EmbryoGrowthTRN=data$EmbryoGrowthTRN, 
                 embryo.stages=data$embryo.stages, 
                 TSP.borders=data$TSP.borders, 
                 TSP.begin=data$TSP.begin, 
                 TSP.end=data$TSP.end, 
                 tsd=data$tsd, 
                 zero=data$zero, 
                Sexed=data$Sexed, Males=data$Males, 
                Temperatures=data$Temperatures, 
                parallel=data$parallel))

}
