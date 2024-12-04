#'
#' @title Expand dataframe for expressions and "all" dimension levels
#' @description Function to expand dataframe for expressions and "all" dimension levels.
#' @param dfr - dataframe with expressions and dimensions to expand from "all" to levels
#' @param lstAlls - list returned by [alls_GetLevels()], or NULL
#' @param verbose - flag (T/F) to print diagnostic info
#' @return dataframe expanded by lstAlls, then by evaluating expressions
#' @details Expressions are evaluated by column row-wise, by column index. Consequently,
#' columns evaluated later can refer to values in columns previously evaluated.
#' Values in expressions will be coerced to the the column type of the original
#' dataframe.
#' @importFrom rlang data_sym
#' @importFrom stringr fixed
#' @importFrom stringr str_detect
#' @importFrom tidyr unnest_longer
#' @import dplyr
#' @export
expandDataframe<-function(dfr,lstAlls=NULL,verbose=FALSE){
  dfrp = dfr;
  if (!is.null(lstAlls)) {
    #--expand "all"s in dimensions, as necessary
    dfrp = alls_ExpandInDataframe(dfrp,lstAlls);
    if (verbose) message("expanded 'alls'")
  }
  for (col in colnames(dfrp)){
    if (verbose) message("processing column '",col,"'");
    col1 = col;
    col2 = rlang::data_sym(col);
    if (is.character(dfrp[[col1]])&&any(stringr::str_detect(dfrp[[col1]],stringr::fixed("(")))){
      if (verbose) message("\tevaluating column '",col,"'");
      dfrp = dfrp |> dplyr::rowwise() |>
                  dplyr::mutate("{col1}":=ifelse(stringr::str_detect(!!col2,fixed("(")),
                                list(as.character(eval(parse(text=!!col2)))),
                                !!col2)) |>
                  dplyr::ungroup();
      if (verbose) message("\tunnesting column '",col,"'");
      dfrp = dfrp |> tidyr::unnest_longer(all_of(col1));
      if (verbose) message("\tsuccessfully expanded column '",col,"'");
      if (verbose) print(dfrp);
    } else {
      if (verbose) message("\tno expansion necessary")
    }
  }
  return(dfrp);
}

#' @title Convert dataframe columns to numeric
#' @description Function to convert dataframe columns to numeric.
#' @param dfr - dataframe to convert
#' @param cols - names of columns to convert to numeric
#' @param verbose - flag (T/F) to print diagnostic info
#' @return dataframe with requisite columns converted to numeric
#' @details None.
#' @export
convertColsToNum<-function(dfr,cols){
  dfrp = dfr;
  for (col in cols){
    col1 = col;
    col2 = rlang::data_sym(col);
    dfrp = dfrp |> dplyr::mutate("{col1}":=as.numeric(!!col2));
  }
  return(dfrp);
}
