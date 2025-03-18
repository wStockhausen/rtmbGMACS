#--extract function-type parameter info
#'
#' @title Extract function-type parameter info from a list
#' @description Function to extract function-type parameter info from a list.
#' @param lst - function-type parameter info list from a `readParamInfo_`
#' @param dms - appropriate `dms...` DimsMap from `setupModelDims`
#' @param process_type - string describing process (for verbose output)
#' @param xtra_cols - character vector of "extra" columns in the "function" section
#' @param verbose - flag (TRUE/FALSE) to print diagnostic info
#' @return a list (see details)
#' @details The output list has elements
#'
#' \itemize{
#' \item{option - format option}
#' \item{Fcns - list with function-level info}
#' \item{dfrCmbs - dataframe with all parameter combinations in vertical format}
#' \item{dfrUniqCmbs - dataframe with unique parameter combinations in vertical format}
#'}
#'
#' @import dplyr
#'
#' @export
#'
extractInfoType1<-function(lst,
                            dms=NULL,
                            process_type="",
                            xtra_cols=character(0),
                            verbose=TRUE){
  if (verbose) message("starting extractInfoType1 for ",process_type);
  #--create list of all dimension levels in `dms` to convert "all"'s to pop levels----
  ##--listAlls is a list with all individual dimension levels, by individual dimension name
  lstAlls = NULL;
  if (!is.null(dms)) lstAlls = alls_GetLevels(dms,verbose=verbose);

  #--create output list----
  out = list();

    ####--functions----
    if (verbose) message("in extractInfoType1 for ",process_type,": processing functions.")
    fcn_cols = c("y","s","r","x","m","p","z","fcn","fcn_idx");
    if (length(xtra_cols)) fcn_cols = c(fcn_cols,xtra_cols);
    dfrFcns = lst$Fcns$dfr |> expandDataframe(lstAlls,verbose=verbose) |>
                dplyr::select(tidyselect::all_of(fcn_cols));
    #browser();

    ###--create full indexing dataframe----
    dfrCmbs = dfrFcns;

    dfrUniqCmbs = dfrCmbs |>
                    dplyr::select(!c(y, s, r, x, m, p, z)) |>
                    dplyr::distinct();

    ###--create output list----
    out = list(Fcns=list(dfrIdxs=lst$Fcns$dfr,
                         dfrRefLvls=dfrRefLvlsFcns,
                         dfrDims2Idxs=dfrFcns),
               dfrCmbs=dfrCmbs,
               dfrUniqCmbs=dfrUniqCmbs);
  return(out);
}

