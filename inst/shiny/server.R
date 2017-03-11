library(shiny)
package.embryogrowth <- require('embryogrowth')

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  # output$resultsInfo <- renderPrint({print(input$MonthRef)})
  
  output$RMUControls <- renderUI({
    df <- subset(embryogrowth::DatabaseTSD, Species==input$Species)$RMU
    df <- c("All", levels(as.factor(as.character(df[df != ""]))))
    selectInput("RMU", "RMU, region or subspecies"
                , choices=as.list(df)
                , selected = "All", multiple = TRUE,
                selectize = TRUE, width = NULL, size = NULL)
  })
  
  output$distPlot <- renderPlot({
    
    if (!package.embryogrowth) {
      par(mar=c(0, 0, 0, 0))
      plot(x=c(0, 1), y=c(0, 1), axes=FALSE,
           xaxt="n", yaxt="n", main="",
           xlab = "", ylab = "",
           xaxs="i", yaxs="i", type="n")
      text(x = 0.5, y=0.5, labels = "The embryogrowth package is not installed!",
           col="red", cex = 2)
      text(x = 0.5, y=0.3, labels = "Contact your site administrator",
           col="green", cex = 2)
    }
    
    buttonP <- input$goButton
    buttonG <- input$predictButton
    if ((buttonP !=0) | (buttonG != 0)) {
      
      isolate({
        zRMU <- input$RMU
        if (is.null(zRMU)) zRMU <- "All"
        
        if (all(zRMU != "All")) {
          zdf <- subset(DatabaseTSD, (RMU %in% zRMU) & 
                          (Species==input$Species) & (Sexed!=0) & 
                          (!is.na(Sexed)) & ((Note != "Sinusoidal pattern")  | (is.na(Note))))
        } else {
          zdf <- subset(DatabaseTSD, (Species==input$Species) & (Sexed!=0) & (!is.na(Sexed)) 
                        & ((Note != "Sinusoidal pattern") | (is.na(Note))))
        }
        
        zmale <- ifelse(input$Male=="1", TRUE, FALSE)
        iP <- input$P
        iS <- input$S
        
        
        if (input$Temperature == "1") {
          Pinit <- c(P = ifelse(is.na(iP), 29, iP), S = ifelse(is.na(iS), -2, iS), K = 0, K1 = 1, K2 = 0)
          zdf <- subset(zdf, select=c("Incubation.temperature", "Males", "Females", "Intersexes", "Reference"))
          zdf$Intersexes <- as.integer(ifelse(is.na(zdf$Intersexes), rep(0, length(zdf$Intersexes)), zdf$Intersexes))
          zdf <- na.omit(zdf)
          output$Data <- renderTable({zdf})
          if (input$Intersexes == "2") zdf$Males <- zdf$Males + zdf$Intersexes
          nx <- sum((zdf$Males == 0) | (zdf$Females == 0))
          output$resultsInfo <- renderText({nx})
          
          if ((nrow(zdf) != 0) & (nx != nrow(zdf))) {
            ztsd.out <- capture.output(ztsd <- tsd(df=zdf, equation=input$Equation, 
                                                   print=TRUE, males.freq = zmale, 
                                                   parameters.initial = Pinit))
          } else {
            ztsd.out <- ""
            ztsd <- NULL
          }
        } else {
          Pinit <- c(P = ifelse(is.na(iP), 55, iP), S = ifelse(is.na(iS), 2, iS), K = 0, K1 = 1, K2 = 0)
          
          meancal <- is.na(zdf$IP.mean) & !is.na(zdf$IP.min) & !is.na(zdf$IP.max)
          zdf[meancal, "IP.mean"] <- (zdf[meancal, "IP.min"]+zdf[meancal, "IP.max"])/2
          zdf <- subset(zdf, select=c("IP.mean", "Males", "Females", "Intersexes", "Reference"))
          zdf$Intersexes <- as.integer(ifelse(is.na(zdf$Intersexes), rep(0, length(zdf$Intersexes)), zdf$Intersexes))
          zdf <- na.omit(zdf)
          output$Data <- renderTable({zdf})
          if (input$Intersexes == "2") zdf$Males <- zdf$Males + zdf$Intersexes
          nx <- sum((zdf$Males == 0) | (zdf$Females == 0))
          output$resultsInfo <- renderText({nx})
          
          if ((nrow(zdf) != 0) & (nx != nrow(zdf))) {
            ztsd.out <- capture.output(ztsd <- tsd(durations = zdf$IP.mean,
                                                   males = zdf$Males,
                                                   females = zdf$Females,
                                                   equation=input$Equation,
                                                   print=TRUE, males.freq = zmale, 
                                                   parameters.initial = Pinit))
          } else {
            ztsd.out <- ""
            ztsd <- NULL
          }
        }
        
        
        
        if (!is.null(ztsd)) {
          
          refT <- paste0(levels(as.factor(as.character(zdf$Reference))), "\n", collapse = "")
        
        output$references <- renderText({refT})
        
        plot(ztsd, males.freq = zmale)
        
        ztsd.out <- gsub("\\[1\\] \\\"", "", ztsd.out)
        ztsd.out <- gsub("\\\"", "", ztsd.out)
        
        output$resultsInfo <- renderText({paste0(ztsd.out[-1], "\n", collapse = "")})
        } else {
          if (nrow(zdf) == 0) {
          output$references <- renderText({"No available data"})
          output$resultsInfo <- renderText({"No available data"})
          } else {
            output$references <- renderText({"At least one datum with mixed sex ratio must be available"})
            output$resultsInfo <- renderText({"At least one datum with mixed sex ratio must be available"})
          }
          plot(x = 1, y=1, type="n", axes=FALSE, xlab="", ylab="", main="")
        }
      })
      
      isolate({
        if (!is.null(ztsd)) {
          valp <- input$Prediction
          valp <- gsub(",", "#", valp)
          valp <- gsub(" +", "#", valp)
          valp <- gsub("#+", "#", valp)
          valp <- as.numeric(unlist(strsplit(valp, "#")))
          # output$p <- renderPrint({print(dput(ztsd))})
          
          o <- predict(ztsd, temperatures=valp)
          output$Prediction <- renderTable(o)
        } else {
          output$Prediction <- renderTable(data.frame())
        }
      })
    } else {
      plot(x = 1, y=1, type="n", axes=FALSE, xlab="", ylab="", main="")
    }
  })
})