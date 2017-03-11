.fonctiontsdMCMC <- function(..., x) {

# .tsd_fit <- function(par, males, N, temperatures, equation)
  # print(data <- data$data)
  # print(x)
  data <- list(...)
return(getFromNamespace(".tsd_fit", ns="embryogrowth")(par=x, males=data$males, N=data$N, fixed.parameters=data$fixed.parameters, temperatures=data$temperatures, equation=data$equation))

}
