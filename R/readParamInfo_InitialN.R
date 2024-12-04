#--read initial population numbers info
#' @title Read initial N-related parameter info
#' @description Function to read initial N-related parameter info.
#' @param conn - connection (or filename) to read from
#' @param debug - flag to print extra info
#' @return nested list with elements (see **Details**)
#' @details The returned list has elements
#' \itemize{
#'   \item{option - initial N calculation option}
#'   \item{lnTotRec - ln-scale total recruitment (only if `option`="free_parameters_v2")}
#'   \item{dfr - dataframe, with parameter specifications (if `option`!="zero_pop")}
#' }
#'
#' @import dplyr
#' @import purrr
#' @import readr
#' @import stringr
#' @md
#' @export
#'
readParamInfo_InitialN<-function(conn,debug=FALSE){
  lns = purrr::keep(stringr::str_trim(readLines(conn,skipNul=TRUE)),
                    \(x)stringr::str_length(x)>0) |>
          extractLines(start="INITIAL_N",end="END");
  iln = 1;
  #--initial N calculation option----
  ##--options are:
  ##--"zero_pop" - no parameters
  ##--"steady-state_unfished"
  ##--"steady-state_fished"
  ##--"free_parameters_v1" - freely estimated parameters
  ##--"free_parameters_v2" - freely estimated parameters, normalized to
  if (debug) cat("reading option\n");
  lst  = extractTextSection(lns,1,1);
  opt  = stringr::str_trim(stringr::str_split_1(lst$txt,"#"))[1];#--calculation option

  #--parse parameter information (as necessary)----
  out = list(option=tolower(opt));
  if (debug) cat(paste0("Initial N option '",opt,"' selected.\n"));
  if (tolower(opt)=="zero_pop"){
    ##--zero_pop----
    #--nothing else to do
  } else if (tolower(opt)=="steady-state_unfished"){
    ##--steady-state_unfished----
    lst   = extractTextSection(lns,1,lst$end+1);
    nRws = as.numeric(stringr::str_split_1(lst$txt,"#")[1]); #--number of lnR0 parameters
    lst   = extractTextSection(lns,nRws+1,lst$end+1);
    out$dfr = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ lnR0 parameters
  } else if (tolower(opt)=="steady-state_fished"){
    ##--steady-state_fished----
    lst   = extractTextSection(lns,1,lst$end+1);
    nRws = as.numeric(stringr::str_split_1(lst$txt,"#")[1]); #--number of lnRini parameters
    lst   = extractTextSection(lns,nRws+1,lst$end+1);
    out$dfr = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ lnRini parameters
  } else if (tolower(opt)=="free_parameters_v1"){
    ##--free_parameters_v1----
    lst   = extractTextSection(lns,1,lst$end+1);
    nRws = as.numeric(stringr::str_split_1(lst$txt,"#")[1]);    #--number of rows defining lnN parameters
    lst   = extractTextSection(lns,nRws+1,lst$end+1);
    out$dfrLnN0 = readr::read_table(I(lst$txt),col_names=TRUE); #--tibble w/ initial lnN0 parameter definitions
  } else if (tolower(opt)=="free_parameters_v2"){
    ##--free_parameters_v2----
    lst  = extractTextSection(lns,1+1,lst$end+1);
    out$dfrLnTotInitPopSize = readr::read_table(I(lst$txt),col_names=TRUE); #--tibble w/ total initial population size
    lst  = extractTextSection(lns,1+1,lst$end+1);
    out$dfrRef = readr::read_table(I(lst$txt),col_names=TRUE);         #--tibble w/ reference class
    lst   = extractTextSection(lns,1,lst$end+1);
    nRws = as.numeric(stringr::str_split_1(lst$txt,"#")[1]);           #--number of rows defining lnN0 parameters
    lst  = extractTextSection(lns,nRws+1,lst$end+1);
    out$dfrLnN0 = readr::read_table(I(lst$txt),col_names=TRUE);        #--tibble w/ lnN0 parameter definitions
  } else {
    stop(paste0("Initial N option '",opt,"' not recognized."));
  }
  return(out);
}
# dirPrj = rstudioapi::getActiveProject();
# source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
# conn=file.path(dirPrj,"testing/InitialN/inputSpecs_InitialN.txt");
# res = readParamInfo_InitialN(conn,TRUE);


