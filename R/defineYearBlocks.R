#' @title Define year blocks
#' @description Function to create a dataframe defining a set of year blocks.
#' @param list - named list of year blocks
#' @param dimsMdl - optional model dimensions tibble
#' @return dataframe with columns "yblk","y"
#' @details Creates a dataframe defining a set of year blocks. If `dimsMdl` is
#' not null, `dimsMdl` is inner-joined to this datafrme, mapping the year blocks to
#' the model dimensions.
#' @examples
#' tbAll = list(all=1949:2022) |> defineYearBlocks();
#' tbRec = list(spinup=1949:1974,normal=1975:2022) |> defineYearBlocks();
#' tbNM = list(normal=c(1949:1979,1985:2022),
#'            heightened=1980:1984) |> defineYearBlocks();
#' tbTCFRetM<-list(early=c(1949:1983,1987:1989),
#'              mid  =c(1990:1996),
#'              PR   =c(2005:2009,2013:2015,2017:2018,2020:2022)) |> defineYearBlocks();
#' @importFrom dplyr bind_rows
#' @importFrom dplyr inner_join
#' @importFrom tibble tibble
#' @export
#'
defineYearBlocks<-function(lst,dimsMdl=NULL){
  lstp = list();
  for (nm in names(lst)){
    lstp[[nm]] = tibble::tibble(yblk=nm,y=lst[[nm]]);
    if (!is.null(dimsMdl))
      lstp[[nm]] = dplyr::inner_join(dimsMdl,lstp[[nm]],by="y");
  }
  dfr=dplyr::bind_rows(lstp);
  return(dfr);
}
