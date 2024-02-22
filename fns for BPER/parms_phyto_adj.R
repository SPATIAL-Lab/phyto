#' Revises prior distributions and parameters
#'
#' This function allows users to revise parameters in the default system file which define prior distributions. This is a highly 
#' recommended step as prior distributions should reflect prior knowledge for the specific study interval. Defaults values for 
#' these revisable parameters reflect late Paleocene to early Eocene in mid to low latitude. 
#' 
#' @details The following is a list of the revisable parameters (default values in parentheses) which define the prior distributions: \cr
#' \cr
#' --- S, T, P --- \cr
#' **'tempC.m'** (30) and **'tempC.sd'** (5) are the mean and stdev for temperature in degrees C  \cr
#' **'sal.m'** (35) and **'sal.sd'** (0.5) are the mean and stdev for salinity  \cr
#' **'press.m'** (6) and **'press.sd'** (1) are the mean and stdev for presure in bar \cr
#' \cr
#' --- carbonate chemistry --- \cr
#' **'pco2.m'** (7.75) and **'pH.l'** (7.45) are the upper and lower bounds on uniform distribution for time step 1 in total scale \cr
#' **'carbchem2.m'** (2000) and **'carbchem2.sd'** (150) are the mean and stdev for 2nd carbonate chemistry variable for time step 1*.  \cr
#' \cr
#' Carbonate chemistry variables in the following units: \cr
#' pH = 'total scale' equivalent  \cr
#' DIC = μmol/kg \cr
#' ALK = μmol/kg \cr
#' CO3 = μmol/kg \cr
#' HCO3 = μmol/kg \cr
#' \cr
#' *Note that variable type, 'cc2ndparm.vt', should be specified in arguments in 'phyto_priors' function. These values will be used if
#' prior type, i.e. 'cc2ndparm.pt', is set to 't1' in 'phyto_priors' argument. These values will not be used if prior type,
#' i.e. 'cc2ndparm.pt', is set to 'ts' in 'phyto_priors' argument.
#' @md
#'
#' @param parms2change Provide a vector of character strings specifying which parameters to revise. The revisable parameters are 
#' defined in 'Details'. 
#' 
#' @param change_values Provide a vector of values which correspond to 'priors2change' vector of parameter names. 
#'
#' @returns Returns list 'parms_phyto_adj' which contains all of the values, revised or otherwise, for revisable parameters.
#'
#' @examples
#' parms_phyto_adj(parms2change = parms2change, change_values = change_values)
#'
#' @export
parms_phyto_adj <- function(parms2change = parms2change, 
                            change_values = change_values){
  
  if(length(parms2change) != length(change_values)){
    stop("Length of 'parms2change' vector much match that of 'change_values' vector")
  }
  
  parms_phyto <- system.file("extdata", "parms_in_phyto.R", package = "BPER")
  source(parms_phyto)
  
  for(i in seq_along(parms2change)){
      assign(noquote(parms2change[i]), change_values[i])
  }
  
  parms_phyto_adj <- list("tempC.m" = tempC.m,
                          "tempC.p" = tempC.p,
                          "sal.m" = sal.m,
                          "sal.p" = sal.p,
                          "pco2.m" = pco2.m,
                          "pco2.p" = pco2.p,
                          "d13C.co2.m" = d13C.co2.m,
                          "d13C.co2.p" = d13C.co2.p,
                          "po4.m" = po4.m,
                          "po4.p" = po4.p,
                          "rm.m" = rm.m,
                          "rm.p" = rm.p) 
  
  class(parms_phyto_adj) = "parms_phyto_adj"
  return(parms_phyto_adj)
}
