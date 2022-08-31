#' Annual average of global temperatures simulated by different CMIP5 GCMs at
#' the planetary scale for the period 1971-2099
#'
#' @name predGlobTemp
#' @docType data
#' @usage data(predGlobTemp)
#' @format matrix 129 years x 20 scenarios
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"predGlobTemp"

#' Equally spaced vector of simulated global temperatures over the period 
#' 1971-2099 for the RCP8.5
#'
#' @name predGlobTempUnique
#' @docType data
#' @usage data(predGlobTempUnique)
#' @format vector of length 13
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"predGlobTempUnique"

#' Years 1971-2099 repeated for the 20 scenarios
#'
#' @name predTime
#' @docType data
#' @usage data(predTime)
#' @format matrix 129 years x 20 scenarios
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"predTime"

#' Equally spaced vector of years over the period 1971-2099
#'
#' @name predTimeUnique
#' @docType data
#' @usage data(predTimeUnique)
#' @format vector of length 13
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"predTimeUnique"

#' scenAvail gives the GCM and RCM which have been used for the 20 climate
#' projections
#'
#' @name scenAvail
#' @docType data
#' @usage data(scenAvail)
#' @format data.frame with 20 rows and two columns: GCM and RCM
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"scenAvail"

#' vecYears gives the years corr. to Y, i.e. from 1971 to 2099
#'
#' @name vecYears
#' @docType data
#' @usage data(vecYears)
#' @format vectors of length 129
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
"vecYears"

#' climate projections of mean winter (DJF) temperature over the SREX region CEU
#' simulated by 20 combinations of CMIP5 GCMs and RCMs for the period 1971-2099
#'
#' @name Y
#' @docType data
#' @usage data(Y)
#' @format matrix 129 years x 20 scenarios
#' @author Guillaume Evin \email{guillaume.evin@inrae.fr}
#' @keywords data
#' @references Seneviratne, S. I. et al. Changes in Climate Extremes and their 
#' Impacts on the Natural Physical Environment, in: Managing the Risks of 
#' Extreme Events and Disasters to Advance Climate Change Adaptation: Special 
#' Report of the Intergovernmental Panel on Climate Change, edited by: Field, 
#' C., Barros, V., Stocker, T., and Dahe, Q., Cambridge University Press, 
#' Cambridge, 109-230, https://doi.org/10.1017/CBO9781139177245.006, 2012
"Y"