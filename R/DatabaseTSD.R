#' Database of TSD information for marine turtles
#' @title Database of TSD information for reptiles
#' @author Marc Girondot \email{marc.girondot@@universite-paris-saclay.fr}
#' @docType data
#' @name DatabaseTSD
#' @encoding UTF-8
#' @description Database of TSD information for reptiles\cr
#' The columns are:\cr
#' \itemize{
#'   \item \code{Species}: Name of the species in binominal nommenclature
#'   \item \code{Subspecies}: Name of the subspecies
#'   \item \code{Country}: From which country the eggs come from
#'   \item \code{Area}: Name of the beach or region the eggs come from
#'   \item \code{RMU}: For marine turtles, name of the RMU for this population; see Wallace, B.P., DiMatteo, A.D., Hurley, B.J., Finkbeiner, E.M., Bolten, A.B., Chaloupka, M.Y., Hutchinson, B.J., Abreu-Grobois, F.A., Amorocho, D., Bjorndal, K.A., Bourjea, J., Bowen, B.W., Duenas, R.B., Casale, P., Choudhury, B.C., Costa, A., Dutton, P.H., Fallabrino, A., Girard, A., Girondot, M., Godfrey, M.H., Hamann, M., Lopez-Mendilaharsu, M., Marcovaldi, M.A., Mortimer, J.A., Musick, J.A., Nel, R., Seminoff, J.A., Troeng, S., Witherington, B., Mast, R.B., 2010. Regional management units for marine turtles: a novel framework for prioritizing conservation and research across multiple scales. Plos One 5, e15465.
#'   \item \code{Incubation.temperature}: Nominal incubation temperature
#'   \item \code{Incubation.temperature.Constant}: Does the incubation temperature was set as constant or CTE was reported
#'   \item \code{Incubation.temperature.Accuracy}: What is the accuracy of the measure of temperature
#'   \item \code{Incubation.temperature.SD}: Experimental SD of incubation temperatures
#'   \item \code{Incubation.temperature.Amplitude}: How much the temperature could fluctuate around nominal temperature
#'   \item \code{Correction.factor}: Difference between the incubator temperature and the eggs temperature
#'   \item \code{IP.min}: Shorter incubation period
#'   \item \code{IP.max}: Longer incubation period
#'   \item \code{IP.mean}: Mean incubation periods
#'   \item \code{IP.SD}: Standard deviation for incubation periods
#'   \item \code{Total}: Total number of eggs incubated
#'   \item \code{Hatched}: Number of hatchlings
#'   \item \code{NotHatched}: Number of embryos with development visible but dead during incubation
#'   \item \code{Undeveloped}: Number of embryos showing no development
#'   \item \code{Intersexes}: Number of individuals intersexes or ambiguous for sex phenotype
#'   \item \code{Males}: Number of individuals indentified as males
#'   \item \code{Females}: Number of individuals indentified as females
#'   \item \code{Sexed}: Number of sexed individuals
#'   \item \code{Clutch}: Identity or number of clutches
#'   \item \code{Reference}: Bibliographic reference
#'   \item \code{Note}: Diverse information for this incubation
#'   \item \code{Version}: Date of the last modification for each record
#' }
#' @keywords datasets
#' @family Functions for temperature-dependent sex determination
#' @usage DatabaseTSD
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(DatabaseTSD)
#' DatabaseTSD.version()
#' totalIncubation_Lo <- subset(DatabaseTSD, 
#'          Species=="Lepidochelys olivacea" & (!is.na(Sexed) & Sexed!=0), 
#'          select=c("Males", "Females", "Incubation.temperature"))
#' tot_Lo <- with(totalIncubation_Lo, tsd(males=Males, females=Females, 
#'  temperatures=Incubation.temperature), parameters.initial = c(P=30.5, S=-0.4))
#'  predict(tot_Lo)
#' }
#' @format A dataframe with raw data.
NULL
