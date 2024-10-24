#'
#' @title Read an ADMB gmacs prj file
#' @description Function to read an ADMB gmacs prj file.
#' @param fn - filename for admb gmacs prj file
#' @param nFlts - number of fleets defined
#' @return list of options for projections, with attribute "type" = "admb prj file inputs"
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
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_split_1
#' @export
#'
readADMB_PrjFile<-function(fn,nFlts){
  lst1 = list(); attr(lst1,"type") = "admb prj file inputs";
  if (file.exists(fn)){
    lns = readLines(fn);
    iln=1;
    splt = stringr::str_split_1(lns[iln],"#")
      #--read MSY option
      iln = findNextLine(lns,iln);
      lst1$optMSY = parseVal(stringr::str_split_1(lns[iln],"#")[1]);

      #--read option to keep F fixed for OFL calculation, by fleet
      lst = list();
      for (f in 1:nFlts){
        #--testing: f = 1;
        iln = findNextLine(lns,iln);
        splt = stringr::str_split_1(lns[iln],"#");
        lst[[f]] = list(f=splt[2],opt=parseVal(splt[1]));
        iln = iln+1;
      }
      lst1$dfrOptsFixedF = dplyr::bind_rows(lst); rm(lst);

      #--read time blocks for OFL
      lst=list();
      for (i in 1:7){
        #--testing: i = 1;
        iln = findNextLine(lns,iln);
        splt = stringr::str_split_1(lns[iln],"#");
        lst[[i]] = list(lbl=splt[2],opt=list(parseVal(splt[1])));
        iln = iln+1;
      }
      lst1$dfrOptsTBs = dplyr::bind_rows(lst); rm(lst);

      #--OFL specification options
      lst=list();
      for (i in 1:7){
        #--testing: i = 1;
        iln = findNextLine(lns,iln);
        splt = stringr::str_split_1(lns[iln],"#");
        lst[[i]] = list(lbl=splt[2],opt=list(parseVal(splt[1])));
        iln = iln+1;
      }
      lst1$dfrOptsOFL = dplyr::bind_rows(lst); rm(lst);

      #--Projection options
      lst=list();
      for (i in 1:22){
        #--testing: i = 1;
        iln = findNextLine(lns,iln);
        splt = stringr::str_split_1(lns[iln],"#");
        lst[[i]] = list(lbl=splt[2],opt=list(parseVal(splt[1])));
        iln = iln+1;
      }
      lst1$dfrOptsPrj = dplyr::bind_rows(lst); rm(lst);

      #--State harvest strategy options
      lst=list();
      for (i in 1:5){
        #--testing: i = 1;
        iln = findNextLine(lns,iln);
        splt = stringr::str_split_1(lns[iln],"#");
        lst[[i]] = list(lbl=splt[2],opt=list(parseVal(splt[1])));
        iln = iln+1;
      }
      lst1$dfrOptsHSs = dplyr::bind_rows(lst); rm(lst);

      #--number of mcmc draws
      iln = findNextLine(lns,iln);
      lst1$optMaxMCMC = parseVal(stringr::str_split_1(lns[iln],"#")[1]);
      iln = iln+1;

      #--full diagnostics option
      iln = findNextLine(lns,iln);
      lst1$optFullDiagnostics = parseVal(stringr::str_split_1(lns[iln],"#")[1]);
  } else {
    stop(paste0("File '",fn,"' does not exist!"));
  }
  return(lst1);
}
