#'
#' @title Get maps from "all" dimension levels to individual levels
#' @description Function to get maps from "all" dimension levels to individual levels.
#' @param dims - list with model dimensions
#' @return list with a dataframe for each dimension mapping "all" to every value.
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @importFrom tidyr expand_grid
#' @export
alls_GetLevels<-function(dims,verbose=FALSE){
  if (verbose) message("in alls_GetLevels.");
  lstLvls = attr(dims,"dmlvs",exact=TRUE);
  lst = list();
  for (lvl in names(lstLvls)){
    lst[[lvl]] = dplyr::bind_rows(
                   tidyr::expand_grid(value=lstLvls[[lvl]],level="all"),
                   tibble::tibble(value=lstLvls[[lvl]],level=lstLvls[[lvl]])
                 );
    names(lst[[lvl]])[2] = lvl;
  }
  if(verbose) print(lst);
  return(lst)
}

#'
#' @title Expand "all" levels in a dataframe
#' @description Function to expand "all" levels in a dataframe.
#' @param dfr - dataframe with dimensions to expand from "all" to levels
#' @param lstAlls - list returned by [alls_GetLevels()]
#' @return expanded dataframe
#' @importFrom rlang data_sym
#' @import dplyr
#' @export
alls_ExpandInDataframe<-function(dfr,lstAlls,verbose=FALSE){
  if (verbose) message("in alls_ExpandInDataframe.");
  dfrp = dfr;
  lvls  = names(lstAlls);
  lvlsp = lvls[lvls %in% names(dfr)];
  lvlsa = lvls[!(lvls %in% names(dfr))];
  #--process levels present
  if (verbose) message("\tprocessing levels present in dfr");
  for (lvl in lvlsp){
    if (verbose) message(paste0("\tprocessing level '",lvl,"'."));
    lvl_sym = rlang::data_sym(lvl);
    dfrp = dfrp |> dplyr::mutate("{lvl}":=as.character(!!lvl_sym)) |>
              dplyr::left_join(lstAlls[[lvl]],relationship="many-to-many");
    if (verbose) {
      message("after left-join.")
      print(dfrp);
    }
    dfrp = dfrp |>
              dplyr::mutate("{lvl}":=ifelse(is.na(value),!!lvl_sym,value)) |>
              dplyr::select(!value);
    if (verbose) {
      message("after mutate and dropping 'value' column")
      print(dfrp);
    }
  }
  #--process levels absent from dfr
  for (lvl in lvlsa){
    lvl_sym = rlang::sym(lvl);
    dfrpp = lstAlls[[lvl]] |> dplyr::select(tidyselect::all_of(lvl)) |> dplyr::filter(tolower(!!lvl_sym)!="all");
    dfrp = dfrp |> dplyr::cross_join(dfrpp)
  }
  #--rearrange columns so "dimensions" are last
  cols = c(names(dfr)[!(names(dfr) %in% names(lstAlls))],names(lstAlls));
  dfrp = dfrp |> dplyr::select(tidyselect::all_of(cols));
  if (verbose) {
    message("finished dfr:");
    print(dfrp);
  }
  return(dfrp);
}
#lstAlls = alls_GetLevels(dims);
#alls_ExpandInDataframe(lst$lstPrcs$dfr,lstAlls[c("r","x","m","a","z")]);
