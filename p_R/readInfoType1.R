#--read model info section
#' @title Read a model info section
#' @description Function to read a model info section
#' @param conn - connection (or filename) to read from
#' @param verbose - flag to print extra info
#' @return nested list with elements (see **Details**)
#' @details The returned list has elements
#' \itemize{
#'   \item{Fcns - list with elements}
#'   \itemize{
#'     \item{n - number of functions defined}
#'     \item{dfr - dataframe with function specs}
#'   }
#' }
#'
#' @import dplyr
#' @import purrr
#' @import readr
#' @import stringr
#' @md
#' @export
#'
readInfoType1<-function(lns,verbose=FALSE){
  iln = 1;
  #--parse input format option----
  if (verbose) cat("Starting readInfoType1\n");
  # lst     = extractTextSection(lns,1,1);
  # option  = (stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--input format option
  #--parse text----
    ###--parse function definitions----
    if (verbose) cat("reading function info\n");
    #lst     = extractTextSection(lns,1,lst$end+1);
    lst     = extractTextSection(lns,1,1);
    nFcns   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of functions
    lst     = extractTextSection(lns,nFcns+1,lst$end+1);
    dfrFcns = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ function defs
    if (verbose) print(dfrFcns);

    ###--create output list----
    if (verbose) cat("creating output list\n");
    out = list(Fcns=list(n=nFcns,dfr=dfrFcns)        #--functions
               );
  if (verbose) cat("Finished readInfoType1\n");
  return(out);
}

