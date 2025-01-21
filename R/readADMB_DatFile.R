
#'
#' @title Read an ADMB gmacs.dat file
#' @description Function to read an ADMB gmacs.dat file (*not* an ADMB gmacs model data file).
#' @param fn - filename for the admb gmacs.dat file
#' @param nFlts - number of fleets defined
#' @return list of options for running an ADMB gmacs model, with attribute "type" = "admb dat file inputs"
#' @details The returned list has the following elements:
#' \itemize{
#'   \item {optMSY}
#'   \item {dfrOptsFixedF - dataframe indicating which fleets have fixed Fs for OFL calculations}
#'   \item {dfrOptsTBs - dataframe with time blocks for OFL calculations}
#'   \item {dfrOptsOFL - dataframe with OFL calculation options}
#'   \item {dfrOptsPrj - dataframe with projection options}
#'   \item {dfrOptsHSs - dataframe with state harvest strategy options}
#'   \item {optMaxMCMC - max number of MCMC draws}
#'   \item {optFullDiagnostics - flag to print full diagnostics}
#' }
#'
#' If it doesn't already exist, a counter `iln` is created in the global environment to keep
#' track of parsing progress. If it already exists, a copy of the current value is made and
#' `iln` is set to 1. Upon function exit, `iln` is reset to its original value if it existed before the function was called;
#' otherwise it is removed.
#'
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_split_1
#' @export
#' @md
#'
readADMB_DatFile<-function(fn,nFlts){
  lst1 = list(); attr(lst1,"type") = "admb dat file inputs";
  if (file.exists(fn)){
    lns = readLines(fn);
    old_iln = NA;
    if (exists("iln",where=.GlobalEnv)) old_iln = iln;
    iln<<-1;
    #--get input files
    lst1$`input files` = parseLines(lns,names=c("data","ctl","prj"));
    #--get units
    lst1$units = parseLines(lns,names=c("wgt","num"));
    #--get stock name
    lst1$stock = parseLines(lns,names=c("stock"));
    #--get jitter info
    lst1$`jitter info` = parseLines(lns,names=c("jitter info"));
    #--flags to output variances
    lst1$`variances to export` = parseLines(lns,
                                            names=c("reference points",
                                                    "recruits",
                                                    "SSB",
                                                    "Fbar",
                                                    "OutDynB0"));
    #--retrospective peel
    lst1$`retro peel` = parseLines(lns,names=c("retro peel"));
    #--other quantities
    lst1$`other` = parseLines(lns,names=c("max phase",
                                          "max num fcn calls",
                                          "calc reference points",
                                          "use pin file",
                                          "verbose"));
    if (!is.na(old_iln)) {
      iln<<-old_iln;#--reset iln to original value, if necessary
    } else {
      rm(iln,envir=.GlobalEnv);
    }
  } else {
    stop(paste0("File '",fn,"' does not exist!"));
  }
  return(lst1);
}
