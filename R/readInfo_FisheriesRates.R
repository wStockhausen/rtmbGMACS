#--read fisheries rates info
#' @title Read fisheries rates info
#' @description Function to read fisheries rates info.
#' @param conn - connection (or filename) to read from
#' @param verbose - flag to print extra info
#' @return nested list with elements (see **Details**)
#' @details
#' @examples
#' # example code
#' if (FALSE){
#'   conn="inputSpecs_FisheriesRates.txt";
#'   lstResults = readInfo_FisheriesRates(conn,verbose=TRUE);
#' }
#'
#' @import dplyr
#' @import purrr
#' @import readr
#' @import stringr
#' @md
#' @export
#'
readInfo_FisheriesRates<-function(conn,verbose=FALSE){
  lns = purrr::keep(stringr::str_trim(readLines(conn,skipNul=TRUE)),
                    \(x)stringr::str_length(x)>0) |>
          extractLines(start="FISHERIES_RATES",end="END");
  out = readInfoType1(lns,verbose);
  return(out);
}

if (FALSE){
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/readInfoType1.R"))
  conn=file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FisheriesRates.txt");
  res = readParamInfo_FisheriesRates(conn,TRUE);
  View(res$Fcns);
 }

