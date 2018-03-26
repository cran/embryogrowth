.fonctionSTRNMCMC <- function(data, x) {

# .STRN_fit <- function(par, EmbryoGrowthTRN, tsd, Sexed, Males, Temperatures) {
  
return(.STRN_fit(par=x, 
                 fixed.parameters=data$fixed.parameters, 
                 EmbryoGrowthTRN=data$EmbryoGrowthTRN, tsd=data$tsd, 
                 embryo.stages=data$embryo.stages, 
                Sexed=data$Sexed, Males=data$Males, 
                Temperatures=data$Temperatures, 
                parallel=data$parallel))

}
