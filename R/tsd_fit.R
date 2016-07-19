.tsd_fit <- function(par, fixed.parameters=NULL, males, N, temperatures, equation) {

#  print(dput(par))
  par <- c(par, fixed.parameters)
  p <- getFromNamespace(".modelTSD", ns="embryogrowth")(par, temperatures, equation)
  p <- ifelse(p<=1E-9, 1E-9, p)
  p <- ifelse(p>=1-1E-9, 1-1E-9, p)
  
#  print(p)

    if (any(is.infinite(p)) | any(is.na(p))) {return(Inf)} else {
      pr <- dbinom(males, N, p, log = TRUE)
#      print(pr)
#      print(-sum(pr))
     return(-sum(pr))
  }
  
}
