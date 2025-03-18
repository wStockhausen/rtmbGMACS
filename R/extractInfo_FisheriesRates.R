#--extract fisheries rates information
#'
#' @title Extract fisheries rates information from info list
#' @description Function to extract fisheries rates information from info list.
#' @param lst - parameter info list from [readInfo_FisheriesRates()]
#' @param dims - `dims` list from `setupModelDims`
#' @return a list (see details)
#' @details The output list has elements
#' \itemize{
#' }
#'
#' @import dplyr
#'
#' @export
#'
extractInfo_FisheriesRates<-function(lst,
                                     dims=NULL,
                                     verbose=TRUE){
  if (verbose) message("starting extractInfo_FisheriesRates.")
  if (FALSE){
    #--NOTE: run this section if you are just stepping through the code for development purposes
    ##--assumes you've created `dims` via `dims = setupModelDims();`
    verbose = TRUE;
    lst = res;#--assumed output from `res  = readInfo_FisheriesRates(conn,verbose=FALSE);`
  }

  #--expand Surveys parameter information----
  ##--inputs are info definitions
  out = extractInfoType1(lst,
                         dims$dmsYSC,
                         "FisheriesRates",
                         xtra_cols=c("flt","cap_idx","ret_idx"),
                         verbose=verbose);
  out$flts = lst$flts;
  return(out);
}

if (FALSE){
  #--test----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R/readInfoSectionType1.R"))
  source(file.path(dirPrj,"R/readInfo_FisheriesRates.R"))
  source(file.path(dirPrj,"R/extractInfoType1.R"))
  source(file.path(dirPrj,"R/extractInfo_FisheriesRates.R"))
  source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))
  dims = setupModelDims();
  conn = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FisheriesRates.txt");
  res1  = readInfo_FisheriesRates(conn,verbose=FALSE);
  res2 = extractInfo_FisheriesRates(res1,dims,verbose=TRUE);
  View(res2$Fcns);
  View(res2f$dfrUniqCmbs);
  View(res2f$dfrUHCs);
}
