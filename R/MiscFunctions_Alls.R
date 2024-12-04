#'
#' @title Get maps from "all" dimension levels to individual levels
#' @description Function to get maps from "all" dimension levels to individual levels.
#' @param dims - list with model dimensions
#' @return list with a dataframe for each dimension mapping "all" to every value.
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @importFrom tidyr expand_grid
#' @export
alls_GetLevels<-function(dims){
  lstLvls = attr(dims$dmsN,"dmlvs",exact=TRUE);
  lst = list();
  for (lvl in names(lstLvls)){
    lst[[lvl]] = dplyr::bind_rows(
                   tidyr::expand_grid(value=lstLvls[[lvl]],level="all"),
                   tibble::tibble(value=lstLvls[[lvl]],level=lstLvls[[lvl]])
                 );
    names(lst[[lvl]])[2] = lvl;
  }
  return(lst)
}

#'
#' @title Expand "all" levels in a dataframe
#' @description Function to expand "all" levels in a dataframe.
#' @param dfr - dataframe with dimensions to expand from "all" to levels
#' @param lstAlls - list returned by [alls_GetLevels()]
#' @return expanded dataframe
#' @importFrom dplyr
#' @export
alls_ExpandInDataframe<-function(dfr,lstAlls){
  dfrp = dfr;
  for (lvl in names(lstAlls)){
    dfrp = dfrp |>
              dplyr::left_join(lstAlls[[lvl]],relationship="many-to-many") |>
              dplyr::mutate("{lvl}":=value) |>
              dplyr::select(!value);
  }
  return(dfrp);
}
#lstAlls = alls_GetLevels(dims);
#alls_ExpandInDataframe(lst$lstPrcs$dfr,lstAlls[c("r","x","m","a","z")]);
