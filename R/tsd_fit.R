.tsd_fit <- function(par, fixed.parameters=NULL, males, N, temperatures, equation) {
  
  #  print(dput(par))
  par <- c(par, fixed.parameters)
  p <- getFromNamespace(".modelTSD", ns="embryogrowth")(par, temperatures, equation)
  p <- ifelse(p<=1E-9, 1E-9, p)
  p <- ifelse(p>=1-1E-9, 1-1E-9, p)
  
  #  print(p)
  
  if (any(is.infinite(p)) | any(is.na(p))) {return(Inf)} else {
    if (any(names(par) == "n")) {
      # sd <- sqrt((males/N)*(1-(males/N))/par["n"])
      if (par["n"] <= 0) {
        return(Inf)
      } else {
        sd <- sqrt((p)*(1-(p))/par["n"])
        pr <- dnorm(males/N, mean=p, sd=sd, log=TRUE)
      }
    } else {
      pr <- dbinom(males, N, p, log = TRUE)
    }
    #      print(pr)
    #      print(-sum(pr))
    L <- -sum(pr)
    attributes(L) <- list(WAIC=pr)
    return(L)
  }
  
}
