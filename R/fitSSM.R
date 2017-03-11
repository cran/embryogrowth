# return the sun of square of differences between the two models

.fitSSM <- function(par, temperatures, growth.rate, fixed.parameters) {
  growth.rate2 <- getFromNamespace(".SSM", ns="embryogrowth")(273.15+temperatures, c(par, fixed.parameters))[[1]]*1E5
  return(sum((growth.rate-growth.rate2)^2, na.rm = TRUE))
}

