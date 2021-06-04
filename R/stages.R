#' Database of embryonic development and thermosensitive period of development for sex determination
#' @title Database of of embryonic development and thermosensitive period of development for sex determination
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @docType data
#' @name stages
#' @encoding UTF-8
#' @description Database of embryonic development and thermosensitive period of development for sex 
#' determination.
#' @references Pieau, C., Dorizzi, M., 1981. Determination of temperature sensitive 
#' stages for sexual differentiation of the gonads in embryos of the turtle, 
#' Emys orbicularis. Journal of Morphology 170, 373-382.
#' @references Yntema, C.L., Mrosovsky, N., 1982. Critical periods and pivotal temperatures for 
#' sexual differentiation in loggerhead sea turtles. Canadian Journal of 
#' Zoology-Revue Canadienne de Zoologie 60, 1012-1016.
#' @references Kaska, Y., Downie, R., 1999. Embryological development of sea turtles (Chelonia mydas, 
#' Caretta caretta) in the Mediterranean. Zoology in the Middle East 19, 55-69.
#' @references Greenbaum, E., 2002. A standardized series of embryonic stages for the emydid 
#' turtle Trachemys scripta. Canadian Journal of Zoology-Revue Canadienne de 
#' Zoologie 80, 1350-1370.
#' @references Magalhães, M.S., Vogt, R.C., Sebben, A., Dias, L.C., de Oliveira, M.F., 
#' de Moura, C.E.B., 2017. Embryonic development of the Giant South American River Turtle, 
#' Podocnemis expansa (Testudines: Podocnemididae). Zoomorphology.
#' @keywords datasets
#' @family Functions for temperature-dependent sex determination
#' @usage stages
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(stages)
#' names(stages)
#' levels(as.factor(stages$Species))
#' # Version of database
#' stages$Version[1]
#' kaska99.SCL <- subset(stages, subset=(Species == "Caretta caretta"), 
#'          select=c("Stage", "SCL_Mean_mm", "SCL_SD_mm", "Days_Begin", "Days_End"))
#' 
#' kaska99.SCL[kaska99.SCL$Stage==31, "Days_Begin"] <- 51
#' kaska99.SCL[kaska99.SCL$Stage==31, "Days_End"] <- 62
#' kaska99.SCL <- na.omit(kaska99.SCL)
#' kaska99.SCL[which(kaska99.SCL$Stage==31), "Stage"] <- c("31a", "31b", "31c")
#' kaska99.SCL <- cbind(kaska99.SCL, 
#'                      Days_Mean=(kaska99.SCL[, "Days_Begin"]+kaska99.SCL[, "Days_End"])/2)
#' kaska99.SCL <- cbind(kaska99.SCL, 
#'                      Days_SD=(kaska99.SCL[, "Days_End"]-kaska99.SCL[, "Days_Begin"])/4)
#' Gompertz <- function(x, par) {
#'    K <- par["K"]
#'    rT <- par["rT"]
#'    X0 <- par["X0"]
#'    y <- abs(K)*exp(log(abs(X0)/abs(K))*exp(-rT*x))
#'    return(y)
#'  }
#' 
#' ML.Gompertz <- function(x, par) {
#'   par <- abs(par)
#'   y <- Gompertz(x, par)
#'   return(sum(-dnorm(y, mean=kaska99.SCL[, "SCL_Mean_mm"], 
#'                     sd=kaska99.SCL[, "SCL_SD_mm"], log=TRUE)))
#' }
#' 
#' parIni <- structure(c(48.66977358, 0.06178453, 0.38640902), 
#'                    .Names = c("K", "rT", "X0"))
#' 
#' fitsize.SCL <- optim(parIni, ML.Gompertz, x=kaska99.SCL[, "Days_Mean"], hessian = TRUE)
#' 
#' # Estimation of standard error of parameters using Hessian matrix
#' sqrt(diag(solve(fitsize.SCL$hessian)))
#' 
#' # Estimation of standard error of parameters using Bayesian  concept and MCMC
#' pMCMC <- structure(list(Density = c("dunif", "dunif", "dunif"), 
#'                         Prior1 = c(0, 0, 0), Prior2 = c(90, 1, 2), 
#'                         SDProp = c(1, 1, 1), 
#'                         Min = c(0, 0, 0), Max = c(90, 1, 2), 
#'                         Init = fitsize.SCL$par), 
#'                    .Names = c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", "Init"), 
#'                    row.names = c("K", "rT", "X0"), class = "data.frame")
#' 
#' Bayes.Gompertz <- function(data, x) {
#'   x <- abs(x)
#'   y <- Gompertz(data, x)
#'   return(sum(-dnorm(y, mean=kaska99.SCL[, "SCL_Mean_mm"], 
#'                     sd=kaska99.SCL[, "SCL_SD_mm"], log=TRUE)))
#' }
#' 
#' mcmc_run <- MHalgoGen(n.iter=50000, parameters=pMCMC, data=kaska99.SCL[, "Days_Mean"], 
#'                      likelihood=Bayes.Gompertz, n.chains=1, n.adapt=100, thin=1, trace=1, 
#'                       adaptive = TRUE)
#' 
#' plot(mcmc_run, xlim=c(0, 90), parameters="K")
#' plot(mcmc_run, xlim=c(0, 1), parameters="rT")
#' plot(mcmc_run, xlim=c(0, 2), parameters="X0")
#' 
#' 1-rejectionRate(as.mcmc(mcmc_run))
#' 
#' par <- mcmc_run$resultMCMC[[1]]
#' 
#' outsp <- t(apply(par, MARGIN = 1, FUN=function(x) Gompertz(0:70, par=x)))
#' 
#' rangqtiles <- apply(outsp, MARGIN=2, function(x) {quantile(x, probs=c(0.025, 0.5, 0.975))})
#' 
#' par(mar=c(4, 4, 2, 1))
#' plot_errbar(x=kaska99.SCL[, "Days_Mean"], y=kaska99.SCL[, "SCL_Mean_mm"], 
#'             errbar.y = 2*kaska99.SCL[, "SCL_SD_mm"], bty="n", las=1, 
#'             ylim=c(0, 50), xlab="Days", ylab="SCL mm", 
#'            xlim=c(0, 70), x.plus = kaska99.SCL[, "Days_End"], 
#'             x.minus = kaska99.SCL[, "Days_Begin"])
#' 
#' lines(0:70, rangqtiles["2.5%", ], lty=2)
#' lines(0:70, rangqtiles["97.5%", ], lty=2)
#' lines(0:70, rangqtiles["50%", ], lty=3)
#' 
#' text(x=50, y=10, pos=4, labels=paste("K=", format(x = fitsize.SCL$par["K"], digits = 4)))
#' text(x=50, y=12.5, pos=4, 
#'    labels=paste("rK=", format(x = fitsize.SCL$par["K"]/39.33, digits = 4)))
#' text(x=50, y=15, pos=4, labels=paste("X0=", format(x = fitsize.SCL$par["X0"], digits = 4)))
#' title("Univariate normal distribution")
#' 
#' # Using a multivariate normal distribution
#' 
#' library(mvtnorm)
#' 
#'  ML.Gompertz.2D <- function(x, par) {
#'    par <- abs(par)
#'   y <- Gompertz(x, par)
#'   L <- 0
#'   for (i in seq_along(y)) {
#'     sigma <- matrix(c(kaska99.SCL$SCL_SD_mm[i]^2, 0, 0, kaska99.SCL$Days_SD[i]^2), 
#'                     nrow=2, byrow=TRUE, 
#'                     dimnames=list(c("SCL_SD_mm", "Days_SD"), c("SCL_SD_mm", "Days_SD")))
#'     L <- L -dmvnorm(x=c(SCL_SD_mm=kaska99.SCL$SCL_Mean_mm[i], 
#'                     Days_SD=kaska99.SCL$Days_Mean[i]), 
#'                     mean= c(SCL_SD_mm=y[i], Days_SD=kaska99.SCL$Days_Mean[i]), 
#'                             sigma=sigma, log=TRUE)
#'   }
#'   return(L)
#' }
#' 
#' parIni <- structure(c(48.66977358, 0.06178453, 0.38640902), 
#'                     .Names = c("K", "rT", "X0"))
#' 
#' fitsize.SCL.2D <- optim(parIni, ML.Gompertz.2D, x=kaska99.SCL[, "Days_Mean"], hessian = TRUE)
#' 
#' # Estimation of standard error of parameters using Hessian matrix
#' sqrt(diag(solve(fitsize.SCL.2D$hessian)))
#' 
#' # Estimation of standard error of parameters using Bayesian  concept and MCMC
#' Bayes.Gompertz.2D <- function(data, x) {
#'   x <- abs(x)
#'   y <- Gompertz(data, x)
#'   L <- 0
#'   for (i in seq_along(y)) {
#'     sigma <- matrix(c(kaska99.SCL$SCL_SD_mm[i]^2, 0, 0, kaska99.SCL$Days_SD[i]^2), 
#'                     nrow=2, byrow=TRUE, 
#'                     dimnames=list(c("SCL_SD_mm", "Days_SD"), c("SCL_SD_mm", "Days_SD")))
#'     L <- L - dmvnorm(x=c(SCL_SD_mm=kaska99.SCL$SCL_Mean_mm[i], 
#'                          Days_SD=kaska99.SCL$Days_Mean[i]), 
#'                     mean= c(SCL_SD_mm=y[i], Days_SD=kaska99.SCL$Days_Mean[i]), 
#'                     sigma=sigma, log=TRUE)
#'   }
#'   return(L)
#' }
#' 
#' pMCMC <- structure(list(Density = c("dunif", "dunif", "dunif"), 
#'                         Prior1 = c(0, 0, 0), Prior2 = c(90, 1, 2), 
#'                         SDProp = c(1, 1, 1), 
#'                         Min = c(0, 0, 0), Max = c(90, 1, 2), 
#'                         Init = fitsize.SCL.2D$par), 
#'                    .Names = c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", "Init"), 
#'                    row.names = c("K", "rT", "X0"), class = "data.frame")
#' mcmc_run.2D <- MHalgoGen(n.iter=50000, parameters=pMCMC, data=kaska99.SCL[, "Days_Mean"], 
#'                      likelihood=Bayes.Gompertz.2D, n.chains=1, n.adapt=100, thin=1, trace=1, 
#'                       adaptive = TRUE)
#' 
#' plot(mcmc_run.2D, xlim=c(0, 90), parameters="K")
#' plot(mcmc_run.2D, xlim=c(0, 1), parameters="rT")
#' plot(mcmc_run.2D, xlim=c(0, 2), parameters="X0")
#' 
#' 1-rejectionRate(as.mcmc(mcmc_run.2D))
#' 
#' par <- mcmc_run.2D$resultMCMC[[1]]
#' 
#' outsp <- t(apply(par, MARGIN = 1, FUN=function(x) Gompertz(0:70, par=x)))
#' 
#' rangqtiles <- apply(outsp, MARGIN=2, function(x) {quantile(x, probs=c(0.025, 0.5, 0.975))})
#' 
#' par(mar=c(4, 4, 2, 1))
#' plot_errbar(x=kaska99.SCL[, "Days_Mean"], y=kaska99.SCL[, "SCL_Mean_mm"], 
#'             errbar.y = 2*kaska99.SCL[, "SCL_SD_mm"], bty="n", las=1, 
#'             ylim=c(0, 50), xlab="Days", ylab="SCL mm", 
#'            xlim=c(0, 70), x.plus = kaska99.SCL[, "Days_End"], 
#'             x.minus = kaska99.SCL[, "Days_Begin"])
#' 
#' lines(0:70, rangqtiles["2.5%", ], lty=2)
#' lines(0:70, rangqtiles["97.5%", ], lty=2)
#' lines(0:70, rangqtiles["50%", ], lty=3)
#' 
#' text(x=50, y=10, pos=4, 
#'      labels=paste("K=", format(x = fitsize.SCL.2D$par["K"], digits = 4)))
#' text(x=50, y=12.5, pos=4, 
#'      labels=paste("rK=", format(x = fitsize.SCL.2D$par["K"]/39.33, digits = 4)))
#' text(x=50, y=15, pos=4, 
#'      labels=paste("X0=", format(x = fitsize.SCL.2D$par["X0"], digits = 4)))
#' title("Multivariate normal distribution")
#' 
#' }
#' @format A list with dataframes including attributes
NULL
