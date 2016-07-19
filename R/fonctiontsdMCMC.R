.fonctiontsdMCMC <- function(data, x) {

# .tsd_fit <- function(par, males, N, temperatures, equation)
  # print(data <- data$data)
  # print(x)
return(getFromNamespace(".tsd_fit", ns="embryogrowth")(par=x, males=data$data$males, N=data$data$N, fixed.parameters=data$data$fixed.parameters, temperatures=data$data$temperatures, equation=data$data$equation))

}
