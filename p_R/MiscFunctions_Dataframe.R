#'
#' @title Expand dataframe for expressions and "all" dimension levels
#' @description Function to expand dataframe for expressions and "all" dimension levels.
#' @param dfrp - dataframe with expressions and dimensions to expand from "all" to levels
#' @param lstAlls - list returned by [alls_GetLevels()], or NULL
#' @param verbose - flag (T/F) to print diagnostic info
#' @return dataframe expanded first by lstAlls, then by evaluating expressions
#' @details Expressions are evaluated by column row-wise, by column index. Consequently,
#' columns evaluated later can refer to values in columns previously evaluated.
#' Values in expressions will be coerced to the the column type of the original
#' dataframe.
#' @importFrom rlang data_sym
#' @importFrom stringr str_detect
#' @importFrom tidyr unnest_longer
#' @import dplyr
#' @export
expandDataframe<-function(dfrp,lstAlls=NULL,verbose=TRUE){
  if (verbose) message("starting expandDataframe.")

  #--expand columns for functional expressions----
  for (col in colnames(dfrp)){
    if (verbose) message("processing column '",col,"'");
    col1 = col;
    col2 = rlang::data_sym(col);
    if (is.character(dfrp[[col1]])){
      if (verbose) {
        message("\tcolumn '",col1,"' is character-valued.");
        print(dfrp[[col1]]);
      }
      if (any(stringr::str_detect(dfrp[[col1]],"\\(|:"))){
        if (verbose) message("\tevaluating column '",col,"'");
        #--NOTE: Need to use `base::ifelse` here:
        #--RTMB defines an `ifelse` that overrides `base::ifelse`, but it throws an error here.
        dfrp = dfrp |> dplyr::rowwise() |>
                    dplyr::mutate("{col1}":=base::ifelse(stringr::str_detect(!!col2,"\\(|:"),
                                                                list(as.character(eval(parse(text=!!col2)))),
                                                                list(as.character(!!col2)))) |>
                    dplyr::ungroup();
        if (verbose) message("\tunnesting column '",col,"'");
        dfrp = dfrp |> tidyr::unnest_longer(all_of(col1));
        if (verbose) {
          message("\tsuccessfully expanded column '",col,"'");
          print(dfrp);
        }
      }
    } else {
      if (verbose) message("\tno expansion necessary")
    }
  }

  #--expand "all"s in dimensions, as necessary----
  if (!is.null(lstAlls)) {
    dfrp = alls_ExpandInDataframe(dfrp,lstAlls,verbose=verbose);
    if (verbose) {
      message("in expandDataframe: expanded 'alls'");
      print(dfrp);
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

#' @title Get the classes (types) of columns in a dataframe
#' @description Function to get the classes (types) of columns in a dataframe.
#' @return a named character vector with values indicating the column class and names indicating the column
#' @export
getColumnTypes<-function(dfr){
  types = sapply(dfr, class);
  return(types);
}

#' @title Get the names of the columns with a given type from a dataframe
#' @description Function to get the names of the columns with a given type from a dataframe.
#' @param dfr - the dataframe
#' @type - the class (as character) of the columns to identify
#' @return character vector of names for columns with class matching `type`
#' @export
getColumnsByType<-function(dfr,type="all"){
  types = sapply(dfr, class);
  if (tolower(type)!="all"){
    types = types[types %in% tolower(type)];
  }
  return(names(types));
}



