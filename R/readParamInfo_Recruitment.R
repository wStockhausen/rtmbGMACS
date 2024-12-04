#--read recruitment info
#' @title Read recruitment-related parameter info
#' @description Function to read recruitment-related parameter info.
#' @param conn - connection (or filename) to read from
#' @param debug - flag to print extra info
#' @return nested list with elements (see **Details**)
#' @details The returned list has elements
#' \itemize{
#'   \item{lstPrcs - process specification list, has elements}
#'   \itemize{
#'     \item{n - number of specifications}
#'     \item{dfr - dataframe with process info}
#'   }
#'   \item{lstTBs - process specification list, has elements}
#'   \itemize{
#'     \item{n - number of time blocks}
#'     \item{dfr - dataframe with time blocks info}
#'   }
#'   \item{lstTotRec - total recruitmnet specification list, has elements}
#'   \item{lstXR - recruitment sex ratio list, has elements}
#'   \item{lstRecZD - recruitment size distribution list, has elements}
#' }
#' where `lstTotRec`, `lstXR`, and `lstRecZD` are each lists having elements
#'   \itemize{
#'     \item{nFcns - number of functions specified}
#'     \item{dfrFcns - dataframe with function info}
#'     \item{nPars - number of main parameters specified}
#'     \item{dfrPars - dataframe with main parameters info}
#'     \item{nRWs - number of random walk parameters specified}
#'     \item{dfrRWs - dataframe with random walk parameters info}
#'     \item{nECs - number of environmental covariates specified}
#'     \item{dfrECs - dataframe with environmental covariates info}
#'   }
#'
#' @import dplyr
#' @import purrr
#' @import readr
#' @import stringr
#' @md
#' @export
#'
readParamInfo_Recruitment<-function(conn,debug=FALSE){
  lns = purrr::keep(stringr::str_trim(readLines(conn,skipNul=TRUE)),
                    \(x)stringr::str_length(x)>0) |>
          extractLines(start="RECRUITMENT",end="END");
  iln = 1;
  #--time blocks----
  if (debug) cat("reading time blocks\n");
  lst  = extractTextSection(lns,1,1);
  nTBs = as.numeric(stringr::str_split_1(lst$txt,"#")[1]);#--number of time blocks
  lst  = extractTextSection(lns,nTBs+1,lst$end+1);
  dfrTBs  = readr::read_table(I(lst$txt),col_names=TRUE);    #--tibble w/ rec time blocks
  lstTBs = list(nTBs=nTBs,dfrTBs=dfrTBs);
  if (debug) print(lstTBs);

  #--total recruitment----
  if (debug) cat("reading total rec\n");
  ##--number of functions----
  if (debug) cat("reading functions\n");
  lst  = extractTextSection(lns,1,lst$end+1);
  nFcns = as.numeric(stringr::str_split_1(lst$txt,"#")[1]); #--number of time series
  ##--functions----
  lst  = extractTextSection(lns,nFcns+1,lst$end+1);
  dfrFcns = readr::read_table(I(lst$txt),col_names=TRUE) |>   #--tibble w/ total recruitment functions
              dplyr::mutate(main_pars=purrr::map(main_pars,\(x)eval(parse(text=x))));
  ##--main parameters----
  if (debug) cat("reading main parameters\n");
  nPars = length(unique(purrr:::list_c(dfrFcns$main_pars))); #--number of main parameters defined
  lst  = extractTextSection(lns,nPars+1,lst$end+1);
  dfrPars = readr::read_table(lst$txt,col_names=TRUE);       #--tibble with main parameter specifications
  ##--RW parameters----
  if (debug) cat("reading RW parameters\n");
  nRWs = sum(parseToDigits(dfrPars$RW_Type),na.rm=TRUE);
  dfrRWs = NULL;
  if (nRWs>0){
    lst  = extractTextSection(lns,nRWs+1,lst$end+1);
    dfrRWs = readr::read_table(lst$txt,col_names=TRUE);       #--tibble with RW parameter specifications
  }
  ##--environmental covariates----
  if (debug) cat("reading env covs\n");
  nECs = sum(dfrPars$Num_Env_Covars);
  dfrECs = NULL;
  if (nECs>0){
    lst  = extractTextSection(lns,nECs+1,lst$end+1);
    stop("TODO: implement environmental covariates section for total recruitment.")
  }
  lstTotRec = list(nFcns=nFcns,dfrFcns=dfrFcns,
                   nPars=nPars,dfrPars=dfrPars,
                   nRWs=nRWs,dfrRWs=dfrRWs,
                   nECs=nECs,dfrECs=dfrECs);
  rm(nFcns,dfrFcns,nPars,dfrPars,nRWs,dfrRWs,nECs,dfrECs);

  #--sex ratio----
  if (debug) cat("reading sex ratio\n");
  ##--number of functions----
  lst  = extractTextSection(lns,1,lst$end+1);
  nFcns = as.numeric(stringr::str_split_1(lst$txt,"#")[1]); #--number of time series
  ##--functions----
  lst  = extractTextSection(lns,nFcns+1,lst$end+1);
  dfrFcns = readr::read_table(lst$txt,col_names=TRUE) |>   #--tibble w/ total recruitment functions
              dplyr::mutate(main_pars=purrr::map(main_pars,\(x)eval(parse(text=x))));
  ##--main parameters----
  nPars = length(unique(purrr:::list_c(dfrFcns$main_pars))); #--number of main parameters defined
  lst  = extractTextSection(lns,nPars+1,lst$end+1);
  dfrPars = readr::read_table(lst$txt,col_names=TRUE);       #--tibble with main parameter specifications
  ##--RW parameters----
  nRWs = sum(parseToDigits(dfrPars$RW_Type),na.rm=TRUE);
  dfrRWs = NULL;
  if (nRWs>0){
    lst  = extractTextSection(lns,nRWs+1,lst$end+1);
    dfrRWs = readr::read_table(lst$txt,col_names=TRUE);       #--tibble with RW parameter specifications
  }
  ##--environmental covariates----
  nECs = sum(dfrPars$Num_Env_Covars);
  dfrECs = NULL;
  if (nECs>0){
    lst  = extractTextSection(lns,nECs+1,lst$end+1);
    stop("TODO: implement environmental covariates section for total recruitment.")
  }
  lstXR = list(nFcns=nFcns,dfrFcns=dfrFcns,
                   nPars=nPars,dfrPars=dfrPars,
                   nRWs=nRWs,dfrRWs=dfrRWs,
                   nECs=nECs,dfrECs=dfrECs);
  rm(nFcns,dfrFcns,nPars,dfrPars,nRWs,dfrRWs,nECs,dfrECs);

  #--recruitment size distribution----
  if (debug) cat("reading rec size dist\n");
  ##--number of functions----
  lst  = extractTextSection(lns,1,lst$end+1);
  nFcns = as.numeric(stringr::str_split_1(lst$txt,"#")[1]); #--number of time series
  ##--functions----
  lst  = extractTextSection(lns,nFcns+1,lst$end+1);
  dfrFcns = readr::read_table(lst$txt,col_names=TRUE) |>   #--tibble w/ total recruitment functions
              dplyr::mutate(main_pars=purrr::map(main_pars,\(x)eval(parse(text=x))));
  ##--main parameters----
  nPars = length(unique(purrr:::list_c(dfrFcns$main_pars))); #--number of main parameters defined
  lst  = extractTextSection(lns,nPars+1,lst$end+1);
  dfrPars = readr::read_table(lst$txt,col_names=TRUE);       #--tibble with main parameter specifications
  ##--RW parameters----
  nRWs = sum(parseToDigits(dfrPars$RW_Type),na.rm=TRUE);
  dfrRWs = NULL;
  if (nRWs>0){
    lst  = extractTextSection(lns,nRWs+1,lst$end+1);
    dfrRWs = readr::read_table(lst$txt,col_names=TRUE);       #--tibble with RW parameter specifications
  }
  ##--environmental covariates----
  nECs = sum(dfrPars$Num_Env_Covars);
  dfrECs = NULL;
  if (nECs>0){
    lst  = extractTextSection(lns,nECs+1,lst$end+1);
    stop("TODO: implement environmental covariates section for total recruitment.")
  }
  lstRecZD = list(nFcns=nFcns,dfrFcns=dfrFcns,
                   nPars=nPars,dfrPars=dfrPars,
                   nRWs=nRWs,dfrRWs=dfrRWs,
                   nECs=nECs,dfrECs=dfrECs);
  rm(nFcns,dfrFcns,nPars,dfrPars,nRWs,dfrRWs,nECs,dfrECs);

  #--process specifications----
  if (debug) cat("reading preocess specs\n");
  lst  = extractTextSection(lns,1,lst$end+1);
  nPrcs = as.numeric(stringr::str_split_1(lst$txt,"#")[1]); #--number of process blocks
  lst   = extractTextSection(lns,nPrcs+1,lst$end+1);
  dfrPrcs = readr::read_table(lst$txt,col_names=TRUE);      #--tibble with process block specifications
  lstPrcs = list(n=nPrcs,dfr=dfrPrcs);
  rm(nPrcs,dfrPrcs);

  #--create return list
  lst = list(lstPrcs=lstPrcs,
             lstTBs=lstTBs,
             lstTotRec=lstTotRec,
             lstXR=lstXR,
             lstRecZD=lstRecZD);
  return(lst);
}

# dirPrj = rstudioapi::getActiveProject();
# source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
# conn=file.path(dirPrj,"testing/Recruitment/inputRecruitmentSpecifications.txt");
# res = readParamInfo_Recruitment(conn,TRUE);
