


.fonctionMCMC <- function(data, x) {

# .fonctionfit<- NULL
# rm(.fonctionfit)


return(.fonctionfit(x, temperatures=data$temperatures, 
                    derivate=data$derivate, weight=data$weight,
                    test=data$test, M0=data$M0, 
                    fixed.parameters=data$fixed.parameters))

}
