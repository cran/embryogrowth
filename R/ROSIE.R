#' Database of TSD information for turtles
#' @title Database of TSD information for turtles
#' @author Caleb Krueger \email{kruegeca@@gmail.com}
#' @docType data
#' @name ROSIE
#' @encoding UTF-8
#' @description Database of incubation information with sex ratio for turtles\cr
#' @keywords datasets
#' @references
#' \insertRef{12889}{embryogrowth}\cr
#' @family Functions for temperature-dependent sex determination
#' @usage ROSIE
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(ROSIE)
#' ROSIE.version()
#' ROSIE <- cbind(ROSIE, RMU.2023=NA)
#' ROSIE[(ROSIE$Species == "Lepidochelys_olivacea") & 
#'       (grepl("Orissa", ROSIE$Location)), 
#'       "RMU.2023"] <- "Northeast Indian"
#' ROSIE[(ROSIE$Species == "Lepidochelys_olivacea") & 
#'       (grepl("Oaxaca", ROSIE$Location) | grepl("Nancite", ROSIE$Location) | 
#'       grepl("Destiladeras", ROSIE$Location) | grepl("Baja California", ROSIE$Location) | 
#'       grepl("Sinaloa", ROSIE$Location) | grepl("Chocó, Colombia", ROSIE$Location) | 
#'       grepl("Jalisco, Mexico", ROSIE$Location)), 
#'       "RMU.2023"] <- "East Pacific"
#' ROSIE[(ROSIE$Species == "Lepidochelys_olivacea") & 
#'       (grepl("Brazil", ROSIE$Location)), 
#'       "RMU.2023"] <- "West Atlantic"
#'       
#' ROSIE[(ROSIE$Species == "Dermochelys_coriacea") & 
#'       (grepl("Suriname", ROSIE$Location) | grepl("French Guiana", ROSIE$Location)), 
#'       "RMU.2023"] <- "Northwest Atlantic"
#' ROSIE[(ROSIE$Species == "Dermochelys_coriacea") & 
#'       (grepl("Playa Grande", ROSIE$Location)), 
#'       "RMU.2023"] <- "East Pacific"
#' ROSIE[(ROSIE$Species == "Dermochelys_coriacea") & 
#'       (grepl("Malaysia", ROSIE$Location)), 
#'       "RMU.2023"] <- "West Pacific"
#'       
#' ROSIE[(ROSIE$Species == "Eretmochelys_imbricata") & 
#'       (grepl("Florida", ROSIE$Location) | grepl("Antigua", ROSIE$Location) | 
#'       grepl("Virgin Islands", ROSIE$Location) | grepl("Puerto Rico", ROSIE$Location) | 
#'       grepl("Campeche", ROSIE$Location) | grepl("Yucatán", ROSIE$Location) | 
#'       grepl("St. Kitts and Nevis", ROSIE$Location)), 
#'       "RMU.2023"] <- "Northwest Atlantic"
#' ROSIE[(ROSIE$Species == "Eretmochelys_imbricata") & 
#'       (grepl("Queensland", ROSIE$Location)), 
#'       "RMU.2023"] <- "Southwest Pacific"
#' ROSIE[(ROSIE$Species == "Eretmochelys_imbricata") & 
#'       (grepl("Brazil", ROSIE$Location)), 
#'       "RMU.2023"] <- "Southwest Atlantic"
#'       
#'       
#' ROSIE[(ROSIE$Species == "Caretta_caretta") & 
#'       (grepl("Georgia", ROSIE$Location) | grepl("South Carolina", ROSIE$Location) | 
#'       grepl("North Carolina", ROSIE$Location) | grepl("Florida", ROSIE$Location) | 
#'       grepl("Texas", ROSIE$Location)), 
#'       "RMU.2023"] <- "Northwest Atlantic"
#' ROSIE[(ROSIE$Species == "Caretta_caretta") & 
#'       (grepl("Cyprus", ROSIE$Location) | grepl("Greece", ROSIE$Location) | 
#'       grepl("Turkey", ROSIE$Location)), 
#'       "RMU.2023"] <- "Mediterranean"
#' ROSIE[(ROSIE$Species == "Caretta_caretta") & 
#'       (grepl("Queensland", ROSIE$Location)), 
#'       "RMU.2023"] <- "South Pacific"
#' ROSIE[(ROSIE$Species == "Caretta_caretta") & 
#'       (grepl("Western Australia", ROSIE$Location)), 
#'       "RMU.2023"] <- "Southeast Indian"
#' ROSIE[(ROSIE$Species == "Caretta_caretta") & 
#'       (grepl("South Africa", ROSIE$Location)), 
#'       "RMU.2023"] <- "Southwest Indian"
#' ROSIE[(ROSIE$Species == "Caretta_caretta") & 
#'       (grepl("Japan", ROSIE$Location)), 
#'       "RMU.2023"] <- "North Pacific"
#' ROSIE[(ROSIE$Species == "Caretta_caretta") & 
#'       (grepl("Brazil", ROSIE$Location)), 
#'       "RMU.2023"] <- "Southwest Atlantic"
#'       
#' ROSIE[(ROSIE$Species == "Chelonia_mydas") & 
#'       (grepl("Queensland", ROSIE$Location)), 
#'       "RMU.2023"] <- "Southwest Pacific"
#' ROSIE[(ROSIE$Species == "Chelonia_mydas") & 
#'       (grepl("Tortuguero", ROSIE$Location) | grepl("British West Indies", ROSIE$Location)), 
#'       "RMU.2023"] <- "North Atlantic"
#' ROSIE[(ROSIE$Species == "Chelonia_mydas") & 
#'       (grepl("Suriname", ROSIE$Location) | grepl("Ascension Island", ROSIE$Location) | 
#'       grepl("Guinea-Bissau", ROSIE$Location)), 
#'       "RMU.2023"] <- "South Atlantic"
#' ROSIE[(ROSIE$Species == "Chelonia_mydas") & 
#'       (grepl("Malaysia", ROSIE$Location) | grepl("Phillipines", ROSIE$Location) | 
#'       grepl("China", ROSIE$Location) | grepl("Taiwan", ROSIE$Location) | 
#'       grepl("Western Australia", ROSIE$Location)), 
#'       "RMU.2023"] <- "East Indian and Southeast Asia"
#' ROSIE[(ROSIE$Species == "Chelonia_mydas") & 
#'       (grepl("Japan", ROSIE$Location)), 
#'       "RMU.2023"] <- "West Central Pacific"
#' ROSIE[(ROSIE$Species == "Chelonia_mydas") & 
#'       (grepl("Cyprus", ROSIE$Location) | grepl("Turkey", ROSIE$Location)), 
#'       "RMU.2023"] <- "Mediterranean"
#' ROSIE[(ROSIE$Species == "Chelonia_mydas") & 
#'       (grepl("Oman", ROSIE$Location)), 
#'       "RMU.2023"] <- "Northwest Indian"
#' ROSIE[(ROSIE$Species == "Chelonia_mydas") & 
#'       (grepl("Hawaii", ROSIE$Location)), 
#'       "RMU.2023"] <- "North Central Pacific"
#'       
#' ROSIE[(ROSIE$Species == "Natator_depressus"), 
#'       "RMU.2023"] <- "nd"
#'       
#' ROSIE[(ROSIE$Species == "Lepidochelys_kempii"), 
#'       "RMU.2023"] <- "Northwest Atlantic"
#'       
#'       
#' totalIncubation_Lo <- subset(ROSIE, 
#'          (Species == "Lepidochelys_olivacea") & (!is.na(Total_Sexed) & Total_Sexed!=0) & 
#'          (Incubation_Setup == "Constant"), 
#'          select=c("Males", "Females", "Mean_Temp", "Latitude", "Longitude", "Location", "RMU.2023"))
#'          
#' library(mapdata)
#' map()
#' scale <- 50
#' title(bquote("Species name: "*italic(.("Lepidochelys olivacea"))))
#' for (l in unique(totalIncubation_Lo$Location)) {
#'   SR_sub <- subset(totalIncubation_Lo, subset = Location == l)
#'   points(x=SR_sub$Longitude[1], y=SR_sub$Latitude[1], pch=19, 
#'   col= 1 + which(l == unique(totalIncubation_Lo$Location)), 
#'   cex=sum(SR_sub$Males + SR_sub$Females, na.rm = TRUE)/scale)
#' }
#' 
#' tot_Lo <- with(totalIncubation_Lo, tsd(males=Males, females=Females, 
#'  temperatures=Mean_Temp), parameters.initial = c(P=30.5, S=-0.4))
#'  plot(tot_Lo, xlim=c(20, 40))
#'  plot(tot_Lo, xlim=c(20, 40), use.ggplot=FALSE)
#'  predict(tot_Lo)
#' }
#' @format A dataframe with raw data.
NULL
