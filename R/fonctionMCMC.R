


.fonctionMCMC <- function(data, x) {

return(info.nests(parameters=x, temperatures=data$temperatures, 
                    derivate=data$derivate, weight=data$weight,
                    test=data$test, M0=data$M0, 
                    fixed.parameters=data$fixed.parameters))

}
