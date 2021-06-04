.STRN_fit <- function(par=NULL                                                                   , 
                      fixed.parameters=NULL                                                      , 
                      equation=NULL                                                              , 
                      tsd=NULL                                                                   , 
                      Sexed=NULL                                                                 , 
                      Males=NULL                                                                 , 
                      sexratio="TSP.TimeWeighted.GrowthRateWeighted.STRNWeighted.sexratio.mean"  , 
                      zero = 1E-9                                                                , 
                      parallel=TRUE                                                              , 
                      NestsResult = NULL                                                         , 
                      embryo.stages=NULL                                                         , 
                      TSP.borders=NULL                                                           , 
                      TSP.begin=NULL                                                             , 
                      TSP.end=NULL                                                               , 
                      out="likelihood"                                                           , 
                      fill=60                                                                    ,
                      verbose = FALSE                                                            ) {
  
  # par=NULL; fixed.parameters=NULL; equation=NULL; TSP.borders=NULL 
  # embryo.stages="Generic.ProportionDevelopment"; TSP.begin=0; TSP.end=0.5
  # EmbryoGrowthTRN=NULL; tsd=NULL; Sexed=NULL; Males=NULL; Temperatures=NULL
  # zero = 1E-9; parallel=TRUE  
  
  # library(embryogrowth)
  # MedIncubation_Cc <- subset(DatabaseTSD, Species=="Caretta caretta" & 
  # RMU=="Mediterranean" & Sexed!=0)
  # tsd <- tsd(males=MedIncubation_Cc$Males, 
  #              females=MedIncubation_Cc$Females, 
  #              temperatures=MedIncubation_Cc$Incubation.temperature, 
  #              par=c(P=29.5, S=-0.1))
  # equation <- NULL 
  # 
  # Males <- c(7, 0, 0, 0, 0, 5, 6, 3, 5, 3, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  # names(Males) <- rev(rev(names(resultNest_4p_SSM$data))[-(1:2)])
  # Sexed <- rep(10, length(males))
  # names(Sexed) <- rev(rev(names(resultNest_4p_SSM$data))[-(1:2)])
  # par <- structure(c(582.567096666926, 2194.0806711639, 3475.28414940385), 
  #                           .Names = c("DHA", "DHH", "T12H"))
  # fixed.parameters <- c(Rho25=100)
  # sexratio <- ""
  # zero <- 1E-9
  # parallel <- TRUE
  # NestsResult <- resultNest_4p_SSM
  # embryo.stages <- "Caretta caretta.SCL"
  # TSP.borders <- NULL  
  # TSP.begin <- 0   
  # TSP.end <- 0.5      
  # sexratio <- "TSP.TimeWeighted.sexratio.mean" 
  
  if (verbose) {
    print(par)
  }
  
  if (!is.null(Sexed)) {
    serafaire <- names(Sexed[(!is.na(Sexed)) & !is.na(Males)])
  } else {
    serafaire <- names(NestsResult$data[1:NestsResult$data[["IndiceT"]]["NbTS"]])
  }
  
  par_fit <- c(par, fixed.parameters)
  
  parSTRN <- par_fit[! (names(par_fit) %in% c("K1", "K2", "K", "P", "S"))]
  
  parTSD <- par_fit[(names(par_fit) %in% c("K1", "K2", "K", "P", "S"))]
  
  if (identical(parTSD, structure(numeric(0), .Names = character(0)))) {
    parTSD <- NULL
  }
  
  if (identical(parSTRN, structure(numeric(0), .Names = character(0)))) {
    parSTRN <- NULL
  }
  
  # Je calcule le SR dépendant de sexratio
  
  if (is.null(tsd)) {
    tsd <- list(par=parTSD, equation=equation)
    class(tsd) <- "tsd"
  }
  if (!is.null(parTSD)) {
    tsd$par <- parTSD
  }
  
  rr <- info.nests(NestsResult=NestsResult                                      ,
                   fixed.parameters=NULL                                        ,
                   series=serafaire                                             ,
                   tsd=tsd                                                      , 
                   embryo.stages=embryo.stages                                  ,
                   TSP.borders=TSP.borders                                      ,
                   TSP.begin=TSP.begin                                          ,
                   TSP.end=TSP.end                                              ,
                   parallel=parallel                                            ,
                   replicate.CI = 0                                             , 
                   out="summary"                                                ,
                   SexualisationTRN = parSTRN                                   ,
                   fill=fill                                                    ,
                   progressbar = FALSE                                          ,
                   warnings = FALSE                                             ,
                   zero=zero                                                    
  )
  sr <- rr$summary[serafaire, sexratio]
  
  if (out=="likelihood") {
    if (is.null(Sexed) | is.null(Males)) {
      return(+Inf)
    } else {
      Sexed <- Sexed[serafaire]
      Males <- Males[serafaire]
      sr <- ifelse(sr == 0, zero, sr)
      sr <- ifelse(sr == 1, 1-zero, sr)
      return(-sum(dbinom(prob=sr, 
                         size=Sexed, x=Males, log=TRUE)))
    }
  } else {
    return(sr)
  }
  
}
