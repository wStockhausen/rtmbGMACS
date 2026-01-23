#' @title Convert CTL file "functions" `tibble` to canonical format
#' @description Convert CTL file "functions" tibble to canonical format.
#' @param tbl - `tibble` to convert
#' @return a tibble with all "functions" ("fcn_idx","fcn_id","function","frame","params","vars","description") columns.
#' @details Converts a "functions" `tibbl`e to canonical format (i.e., al columns). Missing columns
#' in the original tibble are added with `NA`s as values.
#' @import glue
#' @export
#'
getCanonicalFunctionsTbl<-function(tbl){
  cols = c("fcn_idx","fcn_id","function","frame","params","vars","description");
  colsp = cols[!(cols %in% names(tbl))];
  for (col in colsp) tbl = tbl |> tibble::add_column("{col}":=NA_character_);
  return(tbl[,cols]); #--return tbl with columns in canonical order
}
#getCanonicalFunctionsTbl(dfrFcns)

#' @title Convert CTL file "params equations" `tibble` to canonical format
#' @description Convert CTL file "params equations" tibble to canonical format.
#' @param tbl - `tibble` to convert
#' @return a tibble with all "params equations" ("par_idx", "fcn_idx", "par_id", "feEQs", "feCons", "feLink", "reEQs", "reDisp", "reLink", "ecEQs", "ecLink") columns.
#' @details Converts a "params equations" `tibbl`e to canonical format (i.e., al columns). Missing columns
#' in the original tibble are added with `NA`s as values.
#' @import glue
#' @export
#'
getCanonicalParamEquationsTbl<-function(tbl){
  cols = c("par_idx", "fcn_idx", "par_id", "feEQs", "feCons", "feLink", "reEQs", "reDisp", "reLink", "ecEQs", "ecLink");
  colsp = cols[!(cols %in% names(tbl))];
  for (col in colsp) tbl = tbl |> tibble::add_column("{col}":=NA_character_);
  return(tbl[,cols]); #--return tbl with columns in canonical order
}
#getCanonicalParamEquationsTbl(dfrPEQs);
