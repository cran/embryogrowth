#' Result of the mcmc using the nest database using asymmetric normal function
#' @title Result of the fit using the nest database using asymmetric normal function
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name resultNest_mcmc_4p_normal
#' @encoding UTF-8
#' @description Fit using the nest database using asymmetric normal function
#' @references Girondot, M., & Kaska, Y. (2014). A model to predict 
#'             the thermal reaction norm for the embryo growth rate 
#'             from field data. Journal of Thermal Biology, 45, 96-102. 
#'             doi: 10.1016/j.jtherbio.2014.08.005
#' @keywords datasets
#' @usage resultNest_mcmc_4p_normal
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' x <- ChangeSSM(temperatures = (200:350)/10,
#'                parameters = resultNest_4p_SSM4p$par,
#'                initial.parameters = structure(c(3, 7, 11, 32), 
#'                                .Names = c("Scale", "sdL", "sdH", "Peak")), 
#'                control=list(maxit=1000))
#' pfixed <- c(rK=2.093313)
#' resultNest_4p_normal <- searchR(parameters=x$par, fixed.parameters=pfixed, 
#'                          temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#'                          test=c(Mean=39.33, SD=1.92))
#' plotR(resultNest_4p_normal, ylim=c(0, 3))
#' compare_AIC(SSM=resultNest_4p_SSM4p, Asymmetric.normal=resultNest_4p_normal)
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p_normal, accept = TRUE)
#' pMCMC <- structure(list(Density = c("dnorm", "dnorm", "dnorm", "dnorm"
#'            ), Prior1 = c(3, 5, 5, 32), Prior2 = c(5, 2, 2, 2), SDProp = c(0.1, 
#'            0.9, 2, 0.5), Min = c(0, 0, 0, 20), Max = c(100, 20, 20, 40), 
#'            Init = c(2.42695267125171, 5.39959016229435, 4.12171937448323, 
#'            31.746101943379)), .Names = c("Density", "Prior1", "Prior2", 
#'            "SDProp", "Min", "Max", "Init"), row.names = c("Scale", "sdL", 
#'            "sdH", "Peak"), class = "data.frame")
#' resultNest_mcmc_4p_normal <- GRTRN_MHmcmc(result = resultNest_4p_normal, 
#'          n.iter = 50000, parametersMCMC = pMCMC, 
#'          adaptive=TRUE,
#'          intermediate = 1000, filename = "intermediate_mcmc.Rdata")
#' 1-rejectionRate(as.mcmc(resultNest_mcmc_4p_normal))
#' as.parameters(resultNest_mcmc_4p_normal)
#' layout(mat=matrix(1:4, nrow = 2))
#' plot(resultNest_mcmc_4p_normal, parameters = "all", scale.prior = TRUE, las = 1)
#' layout(mat=1)
#' plotR(resultNest_mcmc_4p_normal, ylim=c(0,4), 
#'        main="Asymmetric normal function", show.density=TRUE)
#' }
#' @format A list with fitted information about data(nest)
NULL
