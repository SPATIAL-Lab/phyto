#' Checks and formats phytoplankton proxy data
#'
#' This function allows you to load, QC and format phyto data to use in subsequent 'BPER' functions.
#' Function checks which data are present and whether measurement uncertainties are included.
#' If ages or uncertainties are not included the function stops. Function flags 'incomplete'
#' data object - i.e., missing d13Cpf, or missing Uk'37 - but proceeds. Dependencies: 'readxl'.
#'
#' @param phyto_data User inputs an 11 column R data object and assigns it to 'phyto_data' in arguments. Object columns
#' should be organized as: "site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
#' "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd", "iceco2.data". If 'csv' or 'xlsx' is 
#' to be used instead of R object, provide a character string path to file from directory. 
#' 
#' @param sheet If '.xlsx' is being used provide a character string of the sheet name to read in. 
#'
#' @returns Returns list 'load_phyto'. This list contains 'prox_in' data.frame which has been checked for data completeness
#' and compatibility with other package functions and 'obs_type' which indicates which measurements are included.
#'
#' @examples
#' load_phyto(phyto_data)
#' 
#' @export
load_phyto <- function(phyto_data, 
                       sheet = NULL){
  
  if(is.character(phyto_data)){
    if(grepl(".csv", phyto_data, fixed = TRUE)){
      phyto_data = read.csv(phyto_data)
    } else if(grepl(".xlsx", phyto_data, fixed = TRUE)){
      phyto_data = readxl::read_xlsx(path = file.path, sheet = sheet)
    } else{
      stop("Can only specify a file path to .csv or .xlsx file")
    }
  } else{
    phyto_data = phyto_data
  }
  
  # Check that the correct number of columns are input - i.e. follows the template
  if (ncol(phyto_data) == 11){
    prox_in = data.frame(phyto_data)
    prox_in = prox_in[,c(1:11)]
    names(prox.in) <- c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
                        "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd")
    prox.in <- prox.in[complete.cases(prox.in[,c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
                                                 "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd")]), ]
    prox.in <- transform(prox.in,site.index=as.numeric(factor(site)))
    site.index <- c(prox.in$site.index)
  } else{
    stop("Number of columns in data object do not match template - 11 columns are required even if some are empty")
  }

  # Check that input data frame contains age info, 2 sigma uncertainties are included for each observation
  # and warn if d11B, Mg/Ca or d18O measurements are entirely missing
  if(length(prox_in$age) < 1){
    stop("Must include age data")
  }
  if(length(prox_in$po4.prior) < 1){
    warning("po4.prior values are missing or out of place")
  }
  if(length(prox_in$d13Cmarker.data) < 1){
    warning("d13Cmarker data are missing or out of place")
  }
  if(length(prox_in$d13Cpf.data) < 1){
    warning("d13Cpf data are missing or out of place")
  }
  if(length(prox_in$len.lith.data) < 1){
    warning("len.lith data are missing or out of place")
  }
  if(length(prox_in$Uk.data) < 1){
    warning("Uk data are missing or out of place")
  }
  if(length(prox_in$d13Cmarker.data.sd) < 1 & length(prox_in$d13Cmarker.data) > 1){
    stop("Must include d13Cmarker.data.sd (1 sigma uncertainty) if d13Cmarker data are being used")
  }
  if(length(prox_in$d13Cpf.data.sd) < 1 & length(prox_in$d13Cpf.data) > 1){
    stop("Must include d13Cpf.data.sd (1 sigma uncertainty) if d13Cpf data are being used")
  }
  if(length(prox_in$len.lith.data.sd) < 1 & length(prox_in$len.lith.data) > 1){
    stop("Must include len.lith.data.sd (1 sigma uncertainty) if coccolith length data (len.lith.data) are being used")
  }
  if(length(prox_in$Uk.data.sd) < 1 & length(prox_in$Uk.data) > 1){
    stop("Must include Uk.data.sd (1 sigma uncertainty) if Uk data are being used")
  }

  load_phyto = list("prox_in" = prox_in, "obs_type" = obs_type)
  class(load_phyto) = "load_phyto"
  return(load_phyto)
}
