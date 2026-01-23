#--read selectivity/retention function info
#' @title Read selectivity/retention function info
#' @description Function to read selectivity/retention function info.
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
readParamInfo_SelFcns<-function(conn,debug=FALSE){
  lns = purrr::keep(stringr::str_trim(readLines(conn,skipNul=TRUE)),
                    \(x)stringr::str_length(x)>0) |>
          extractLines(start="SEL_FCNS",end="END");
  iln = 1;
  #--parse function definitions----
  if (debug) cat("reading function info\n");
  lst     = extractTextSection(lns,1,1);
  nFcns   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of functions
  lst     = extractTextSection(lns,nFcns+1,lst$end+1);
  dfrFcns = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ function defs
  lst     = extractTextSection(lns,1,lst$end+1);
  nLvls   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of reference levels
  if (nLvls>0){
    if (debug) cat("reading function reference levels info\n");
    lst    = extractTextSection(lns,nLvls+1,lst$end+1);
    dfrRefLvls = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ function reference levels
  } else {dfrRefLvls=NULL;}

  #--parse parameter information (as necessary)----
  ##--Main parameters----
  if (debug) cat("reading main parameters info\n");
  lst     = extractTextSection(lns,1,lst$end+1);
  nMPs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of main parameters
  if (nMPs>0){
    lst     = extractTextSection(lns,nMPs+1,lst$end+1);
    dfrMPs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ main parameter defs
  } else {dfrMPs=NULL;}

  ##--parse offset parameters----
  if (debug) cat("reading offset parameters info\n");
  lst     = extractTextSection(lns,1,lst$end+1);
  nOPs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of offset parameters
  if (nOPs>0){
    lst     = extractTextSection(lns,nOPs+1,lst$end+1);
    dfrOPs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ offset parameter defs
  } else {dfrOPs=NULL;}

  ##--parse "devs" parameters----
  if (debug) cat("reading devs parameters info\n");
  lst     = extractTextSection(lns,1,lst$end+1);
  nDPs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of devs parameters
  if (nDPs>0){
    lst     = extractTextSection(lns,nDPs+1,lst$end+1);
    dfrDPs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ devs parameter defs
    ###--reference levels----
    lst     = extractTextSection(lns,1,lst$end+1);
    nLvls_DPs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of devs reference level defs
    if (nLvls_DPs>0){
      if (debug) cat("reading devs parameter referece levels info\n");
      lst     = extractTextSection(lns,nLvls_DPs+1,lst$end+1);
      dfrRefLvls_DPs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ devs reference level def
    } else {dfrRefLvls_DPs=NULL;}
  } else {dfrDPs=NULL; nLvls_DPs=0; dfrRefLvls_DPs=NULL;}

  ##--parse environmental covariates----
  if (debug) cat("reading environmental covariates info\n");
  lst     = extractTextSection(lns,1,lst$end+1);
  nECs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of environmental covariates
  if (nECs>0){
    lst     = extractTextSection(lns,nECs+1,lst$end+1);
    dfrECs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ environmental covariate defs
    ###--reference levels----
    lst     = extractTextSection(lns,1,lst$end+1);
    nLvls_ECs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of environmental covariate reference level defs
    if (nLvls_ECs>0){
      if (debug) cat("reading environmental covariate reference levels info\n");
      lst     = extractTextSection(lns,nLvls_ECs+1,lst$end+1);
      dfrRefLvls_ECs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ environmental covariate reference level def
    } else {dfrRefLvls_ECs=NULL;}
  } else {dfrECs=NULL; nLvls_ECs=0; dfrRefLvls_ECs=NULL;}

  #--parse functional priors----
  if (debug) cat("reading functional priors info\n");
  lst     = extractTextSection(lns,1,lst$end+1);
  nFPs   = as.integer(stringr::str_trim(stringr::str_split_1(lst$txt,"#")))[1];#--number of functional priors
  if (nFPs>0){
    lst     = extractTextSection(lns,nFPs+1,lst$end+1);
    dfrFPs = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ environmental covariate defs
  } else {dfrFPs=NULL;}

  #--create output list----
  out = list(Fcns=list(n=nFcns,dfr=dfrFcns,
                       reflvls=list(n=nLvls,dfr=dfrRefLvls)),        #--functions
             MPs=list(n=nMPs,dfr=dfrMPs),                            #--main parameters
             OPs=list(n=nOPs,dfr=dfrOPs),                            #--offset parameters
             DPs=list(n=nDPs,dfr=dfrDPs,
                      reflvls=list(n=nLvls_DPs,dfr=dfrRefLvls_DPs)), #--devs parameters
             ECs=list(n=nECs,dfr=dfrECs,
                      reflvls=list(n=nLvls_ECs,dfr=dfrRefLvls_ECs)), #--env. covars
             FPs=list(n=nFPs,dfr=dfrFPs)                             #--functional priors
             )
  return(out);
}
# dirPrj = rstudioapi::getActiveProject();
# source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
# conn=file.path(dirPrj,"testing/testSelFunctions/inputSpecs_SelFcns.txt");
# res = readParamInfo_SelFcns(conn,TRUE);


