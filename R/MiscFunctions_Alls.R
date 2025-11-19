#'
#' @title Get maps from "all" dimension levels to individual levels
#' @description Function to get maps from "all" dimension levels to individual levels.
#' @param dfrDims - dataframe with model dimensions to use for expansion
#' @return list with a dataframe for each dimension mapping "all" to every value (and every value to itself).
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @importFrom tidyr expand_grid
#' @export
alls_GetLevels<-function(dfrDims,verbose=FALSE){
  if (verbose) message("in alls_GetLevels.");
  lstLvls = attr(dfrDims,"dmlvs",exact=TRUE);
  lst = list();
  for (lvl in names(lstLvls)){
    lst[[lvl]] = dplyr::bind_rows(
                   tidyr::expand_grid(value=lstLvls[[lvl]],level="all"),
                   tibble::tibble(value=lstLvls[[lvl]],level=as.character(lstLvls[[lvl]]))
                 );
    names(lst[[lvl]])[2] = lvl;
  }
  if(verbose) print(lst);
  return(lst)
}

#'
#' @title Expand "all" levels in a dataframe
#' @description Function to expand "all" levels in a dataframe.
#' @param dfr - dataframe with dimensions to expand from "all" to the individual levels
#' @param lstAlls - list in format returned by [alls_GetLevels()]
#' @return expanded dataframe
#' @details The `lstAlls` object can either be a list of dataframes returned from [alls_GetLevels()] mapping
#' "all"s to each model dimension or a "roll your own" list.
#'
#' The "roll your own"
#' list can be in the format returned from [alls_GetLevels()] (i.e., a dataframe of levels for each "dimension" to expand)
#' or as a named list with a vector of the values to expand to for each "dimension". In either case, the names of the list
#' elements correspond to the associated "dimension" to expand.
#'
#' @examplesIf FALSE
#' # example code
#' dfr = tibble::tibble(x=1:3,y="all");
#' #--`lstAlls` in format returned from [alls_GetLevels()]
#' lstAlls = list(y=tibble::tibble(value=as.character(paste0("y",c(1:2,1:2))),y=c("all","all",paste0("y",as.character(1:2)))));
#' alls_ExpandInDataframe(dfr,lstAlls);
#'
#--listAlls in simple vector format
#' lstAlls = list(y=c("y1","y2"));
#' alls_ExpandInDataframe(dfr,lstAlls);
#'
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
    dfrAlls = lstAlls[[lvl]];
    if (is.vector(dfrAlls)) {
      vals = dfrAlls;
      dfrAlls = tibble::tibble(value=c(vals,vals)) |> dplyr::mutate("{lvl}":=c(as.character(vals),rep("all",length(vals))));
    }
    dfrp = dfrp |> dplyr::mutate("{lvl}":=as.character(!!lvl_sym)) |>
              dplyr::left_join(dfrAlls,relationship="many-to-many");
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
    dfrAlls = lstAlls[[lvl]];
    if (is.vector(dfrAlls)) {
      vals = dfrAlls;
      dfrAlls = tibble::tibble(value=c(vals,vals)) |> dplyr::mutate("{lvl}":=c(as.character(vals),rep("all",length(vals))));
    }
    dfrpp = dfrAlls |> dplyr::select(tidyselect::all_of(lvl)) |> dplyr::filter(tolower(!!lvl_sym)!="all");
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
