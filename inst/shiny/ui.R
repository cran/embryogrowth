# library(shiny); library("embryogrowh"); runApp("/Users/marc/Documents/Espace_de_travail_R/shiny/tsd")
# library(shiny); runApp("http://max3.ese.u-psud.fr:3838/phenology/")


library(shiny)
library(embryogrowth)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  
  wellPanel(
    h1("Temperature-dependent sex determination", align = "center")
    , h2("at constant temperature", align = "center")
    , p(em("Temperature-dependent sex determination (TSD) is a type of environmental 
        sex determination in which the temperatures experienced during embryonic 
        development determine the sex of the offspring. It is present in all crocodilians, 
        many turtles and squamates."))
    , p(strong("The embryogrowth R package is a set of tools used 
               to model embryo growth and sexual phenotype linked to 
               temperature."))
    , p("This web server version v. 3.62 is a simplified version of the complete tools available ", 
        a("here."
          , href="https://cran.r-project.org/package=embryogrowth"
          , target="_blank"))
    , p("Embryogrowth package is developped by "
        , a("Marc Girondot"
            , href="http://max2.ese.u-psud.fr/epc/conservation/Girondot/Publications/Marc.html"
            , target="_blank"))
    , p(paste0("Database version ", as.character(DatabaseTSD$Version[1])))
  )
    , wellPanel(
    p(" The methods have been published in:")
    , p(a(
      img(src="PDF.png", height=35, width=35), 
      "Girondot, M. 1999 Statistical description of temperature-dependent 
      sex determination using maximum likelihood. Evolutionary Ecology Research 
      1, 479-486.", 
      href="https://www.researchgate.net/publication/242076781_Statistical_description_of_temperature-dependent_sex_determination_using_maximum_likelihood",
      target="_blank"))
    , p(a(
      img(src="PDF.png", height=35, width=35), 
      "Godfrey, M.H., Delmas, V. & Girondot, M. 2003 Assessment of 
      patterns of temperature-dependent sex determination using maximum likelihood 
      model selection. Ecoscience 10, 265-272.", 
      href="https://www.researchgate.net/profile/Marc_Girondot/publications?sorting=newest&page=5",
      target="_blank"))
    , p(a(
      img(src="PDF.png", height=35, width=35), 
      "Hulin, V., Delmas, V., Girondot, M., Godfrey, M.H. & Guillon, J.-M. 2009 
      Temperature-dependent sex determination and global change: Are some species at 
      greater risk? Oecologia 160, 493-506.", 
      href="https://www.researchgate.net/publication/24192820_Temperature-dependent_sex_determination_and_global_change_are_some_species_at_greater_risk",
      target="_blank"))
    , helpText(strong("These methods must not be used for incubation data at non-constant temperatures, both 
               for temperature versus sex ratio and incubation duration versus sex ratio."))
    )
  
  # Sidebar with a slider input for the number of bins
  , sidebarLayout(
    sidebarPanel(
      wellPanel(
        h4("Load dataset from litterature")
        , selectInput("Species", "Species"
                      , choices=as.list(levels(as.factor(embryogrowth::DatabaseTSD$Species)))
                      , selected = NULL, multiple = FALSE,
                      selectize = TRUE, width = NULL, size = NULL)
        # , selectInput("RMU", "RMU"
        #               , choices=as.list(c("All", levels(as.factor(embryogrowth::DatabaseTSD$RMU))[-1]))
        #               , selected = NULL, multiple = TRUE,
        #               selectize = TRUE, width = NULL, size = NULL)
        , uiOutput("RMUControls")
        , helpText("RMU are the Regional Managment Units defined for marine turles in ", 
                   a("Wallace B.P., et al. (2010) Regional management units for marine turtles: 
                     a novel framework for prioritizing conservation and research across multiple 
                     scales. PLoS One, 5(12), e15465.", 
                     href="https://www.researchgate.net/publication/49772535_Wallace_BP_DiMatteo_AD_Hurley_BJ_Finkbeiner_EM_Bolten_AB_et_al_Regional_management_units_for_marine_turtles_a_novel_framework_for_prioritizing_conservation_and_research_across_multiple_scales_PLoS_ONE",
                     target="_blank"))
        , actionButton("goButton", "Go!")
      )
      
      ,  wellPanel(h4("Options to evaluate TSD")
                   , radioButtons("Male", "Should the data presented as male or female frequency?", list(Male=1, Female=2), selected=1, inline = TRUE)
                   , radioButtons("Temperature", "Analyze sex ratio against temperature or incubation duration ?", list(Temperature=1, Duration=2), selected=1, inline = TRUE)
                   , selectInput("Equation", "What equation should be used?"
                                 , choices=list("GSD", "Logistic", "Hill", "Richards", "Double-Richards", "Hulin")
                                 , selected = "Logistic", multiple = FALSE,
                                 selectize = TRUE, width = NULL, size = NULL)
                   , radioButtons("Intersexes", "How to count intersexes?", list(Discard=1, Males=2), selected=1, inline = TRUE)
                   , h4("Initial values for fit")
                   , numericInput("P", "P value", "", width = "50%")
                   , numericInput("S", "S value", "", width = "50%")
                   , helpText("In most of the cases, you can let these boxes empty.")
      )
      ,  wellPanel(h4("Prediction")
                   , textInput("Prediction", "Values for which sex ratio should be predicted:", "30")
                   , helpText("Values can be separated with space or commas.")
                   , actionButton("predictButton", "Predict")

                   )
    )
    ,
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
      , h4("Synthesis for TSD model")
      , verbatimTextOutput(outputId="resultsInfo")
      , h4("Data used for the fit")
      , tableOutput(outputId="Data")
      # , verbatimTextOutput(outputId="Data")
      , h4("References")
      , verbatimTextOutput(outputId="references")
      , h4("Prediction")
      , tableOutput(outputId="Prediction")
      # , verbatimTextOutput(outputId="Prediction")
      #, verbatimTextOutput(outputId="p")
      
    )
  )
  ))