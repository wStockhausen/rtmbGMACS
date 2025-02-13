#--read natural mortality info
#' @title Read natural mortality info
#' @description Function to read natural mortality info.
#' @param conn - connection (or filename) to read from
#' @param verbose - flag to print extra info
#' @return nested list with elements (see **Details**)
#' @details If `option` = "function", the returned list has elements
#' \itemize{
#'   \item{option}
#'   \item{Fcns - list with elements}
#'   \itemize{
#'     \item{n - number of functions defined}
#'     \item{dfr - dataframe with function specs}
#'     \item{reflvls - list with number of reference levels definied (`n`) and a dataframe with the levels (`dfr`)}
#'   }
#'   \item{MPs - list with main parameters info; elements are}
#'   \itemize{
#'     \item{n - number of main parameters defined}
#'     \item{dfr - dataframe with main parameters specs}
#'   }
#'   \item{OPs - list with offset parameters info; elements are}
#'   \itemize{
#'     \item{n - number of offset parameters defined}
#'     \item{dfr - dataframe with offset parameters specs}
#'   }
#'   \item{DPs - list with devs parameters info; elements are}
#'   \itemize{
#'     \item{n - number of dev parameters defined}
#'     \item{dfr - dataframe with dev parameter specs}
#'     \item{reflvls - list with number of reference levels defined (`n`) and a dataframe with the levels (`dfr`)}
#'   }
#'   \item{REs - list with random effects parameters info; elements are}
#'   \itemize{
#'     \item{n - number of random effects parameters defined}
#'     \item{dfr - dataframe with random effects parameter specs}
#'     \item{reflvls - list with number of reference levels defined (`n`) and a dataframe with the levels (`dfr`)}
#'   }
#'   \item{ECs - list with environmental covariates info; elements are}
#'   \itemize{
#'     \item{n - number of environmental covariates defined}
#'     \item{dfr - dataframe with environmental covariates specs}
#'     \item{reflvls - list with number of reference levels defined (`n`) and a dataframe with the levels (`dfr`)}
#'   }
#'   \item{FPs - list with functional priors info; elements are}
#'   \itemize{
#'     \item{n - number of functional priors defined}
#'     \item{dfr - dataframe with functional priors specs}
#'   }
#' }
#'
#' If `option` = "data", a list with the following elements is returned:
#' \itemize{
#'   \item{option}
#'   \item{transform - name of function to transform input values to kg}
#'   \item{dfr - dataframe with value specs }
#' }
#'
#' @examples
#' # example code reading "vertical" data format
#' if (FALSE){
#'   conn="inputSpecs_NaturalMortality.data_vertical.txt";
#'   lstResults = readParamInfo_NaturalMortality(conn,verbose=TRUE);
#' }
#'
#' # example code reading "horizontal" data format
#' if (FALSE){
#'   conn="inputSpecs_NaturalMortality.data_horizontal.txt";
#'   lstResults = readParamInfo_NaturalMortality(conn,verbose=TRUE);
#' }
#'
#' # example code reading "function" format
#' if (FALSE){
#'   conn="inputSpecs_NaturalMortality.function.txt";
#'   lstResults = readParamInfo_NaturalMortality(conn,verbose=TRUE);
#' }
#'
#' @import dplyr
#' @import purrr
#' @import readr
#' @import stringr
#' @md
#' @export
#'
readParamInfo_NaturalMortality<-function(conn,verbose=FALSE){
  lns = purrr::keep(stringr::str_trim(readLines(conn,skipNul=TRUE)),
                    \(x)stringr::str_length(x)>0) |>
          extractLines(start="NATURAL_MORTALITY",end="END");
  out = readParamInfoSectionType1(lns,verbose);
  return(out);
}

if (FALSE){
  #--test "data"-type vertical input----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/readParamInfoSectionType1.R"))
  conn=file.path(dirPrj,"testing/testNaturalMortality/inputSpecs_NaturalMortality.data-vertical.txt");
  resV = readParamInfo_NaturalMortality(conn,TRUE);
  View(resV$dfr);
}

if (FALSE){
  #--test "data"-type horizontal input----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/readParamInfoSectionType1.R"))
  conn=file.path(dirPrj,"testing/testNaturalMortality/inputSpecs_NaturalMortality.data-horizontal.txt");
  resH = readParamInfo_NaturalMortality(conn,TRUE);
  View(resH$dfr);
}

if (FALSE){
  #--test "function"-type input----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/readParamInfoSectionType1.R"))
  conn=file.path(dirPrj,"testing/testNaturalMortality/inputSpecs_NaturalMortality.function.txt");
  resF = readParamInfo_NaturalMortality(conn,TRUE);
  View(resF);
}

