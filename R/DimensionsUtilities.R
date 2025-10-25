#'
#' @title Convert non-numeric dimension values to factor levels in a dimensions map
#' @description Function to convert non-numeric dimension values to factor levels in a dimensions map.
#' @param dfr - dimensions map dataframe/tibble (e.g. output from [createSparseDimsMap()])
#'
#' @return - dimensions map with non-numeric dimension values converted to factors
#'
#' @details Dimension columns that are already factors are not converted.
#'
#' @export
#'
createDimsFactors<-function(dfr){
    dmnms = attr(dfr,"dmnms");
    dmlvs = attr(dfr,"dmlvs");
    for (dmnm in dmnms){
        if (!(is.factor(dfr[[dmnm]])||is.numeric(dfr[[dmnm]])))
            dfr[[dmnm]] = factor(dfr[[dmnm]],levels=dmlvs[[dmnm]]);
    }
    return(dfr);
}
# dfrp = createDimsFactors(dfrSparse); str(dfrp);
# dfrp = createDimsFactors(dfrp);      str(dfrp);

#'
#' @title Convert values to dimension factor levels in a dataframe
#' @description Function to convert dataframe values to factor levels in a dimensions map.
#' @param dfr - dataframe to convert
#' @param dfrDms - dimensions map dataframe/tibble (e.g. output from [createSparseDimsMap()])
#' @param contrasts - character string, list, FALSE, or NULL defining contrasts to apply to (eventual) factor columns
#' @param all - flag (T/F) to convert all columns to factors (default=TRUE)
#' @return - dataframe with values converted to factors
#'
#' @details If a column name in `dfr` matches a dimension name,
#' factor or character columns are converted to a factor with the same levels as the dimension.
#' This also occurs if the column is numeric and `all` is TRUE. Otherwise, the
#' column is converted to a factor using `unique` to determine the levels
#' if it is a character column or if it numeric and `all` is TRUE.
#'
#' If `contrasts` is `FALSE`, then no contrasts are applied are applied to any factor. If contrasts is NULL,
#' then the default contrast given by options("contrast") is applied  to each factor.
#' Otherwise, contrasts should be a list, with the contrast associated with each name in the list applied to
#' the corresponding factor.
#'
#' @export
#'
convertToDimsFactors<-function(dfr,dfrDms,contrasts=NULL,all=FALSE,verbose=FALSE){
    dmnms = attr(dfrDms,"dmnms");
    dmlvs = attr(dfrDms,"dmlvs");
    # if (verbose) print(dmnms);
    # if (verbose) print(dmlvs);
    # if (verbose) print(contrasts);
    if (!is.null(contrasts)){
      if (is.character(contrasts)){
        if (tolower(contrasts)=="false") {
          if (verbose) cat("setting all contrasts to FALSE\n");
          contrasts = FALSE;
        } else
        if (stringr::str_starts(contrasts,"ident")){
          if (verbose) cat("setting contrasts to NULL so using default:",options(contrasts)$contrasts["unordered"],"\n");
          contrasts = NULL;
        } else {
          if (verbose) cat("evaluating lists of contrasts\n");
          if (stringr::str_starts(contrasts,"list")){
            eval(parse(text=paste0("contrasts=",contrasts)));
          } else {
            eval(parse(text=paste0("contrasts=list(",contrasts,")")));
          }
          if (verbose) cat("class(contrasts) = ",class(contrasts),"\n");
        }
      } #--is.character(contrasts)
    }#--!is.null(contrasts)
    if (verbose) {cat("class(contrasts):",class(contrasts),"\n");}
    if (verbose) {cat("contrasts:\n"); print(contrasts);}
    if (verbose) cat("\napplying contrasts\n")
    for (colnm in names(dfr)){
      if (verbose) cat("checking",colnm,"\n");
      if (colnm %in% dmnms){
        #--column name matches a dimension
        levs = dmlvs[[colnm]];
        if (verbose) cat("\tcolumn name matches a dimension\n");
        if (all){
          # convert column to factor with same levels as dimension
          if (verbose) cat("converting factor to same levels as dimension\n");
          dfr[[colnm]] = factor(as.character(dfr[[colnm]]),levels=levs);
          if (is.list(contrasts)){
            if (colnm %in% names(contrasts)){
              if (!is.logical(contrasts[[colnm]])) {
                dfr[[colnm]] = C(dfr[[colnm]],contrasts[[colnm]]);
              }
            }
          }
        } else {
          #--don't convert numeric columns
          if (!is.numeric(dfr[[colnm]])){
            if (verbose) cat("converting",colnm,"to factor\n")
            dfr[[colnm]] = factor(as.character(dfr[[colnm]]),levels=levs);
            if (verbose) cat(class(dfr[[colnm]]),"\n");
            if (is.list(contrasts)){
              if (colnm %in% names(contrasts)){
                if (!is.logical(contrasts[[colnm]])) {
                  if (verbose) cat("setting",contrasts[[colnm]],"on",colnm,"\n")
                  dfr[[colnm]] = C(dfr[[colnm]],contrasts[[colnm]]);
                }
              }
            }
          }
        }
      } else {
        #--column name does not match a dimension
        if (all&&(!is.factor(dfr[[colnm]]))){
          #--convert column to factor if not already a factor
          levs = unique(as.character(dfr[[colnm]]));
          dfr[[colnm]] = factor(as.character(dfr[[colnm]]),levels=levs);
          if (is.list(contrasts)){
            if (colnm %in% names(contrasts))
              if (!is.logical(contrasts[[colnm]])) dfr[[colnm]] = C(dfr[[colnm]],contrasts[[colnm]]);
          }
        } else {
          if (is.factor(dfr[[colnm]])){
            if (colnm %in% names(contrasts)){
              #--add contrast
              if (verbose) cat("Adding contrast",contrasts[[colnm]],"to",colnm,"\n");
              if (!is.logical(contrasts[[colnm]])) dfr[[colnm]] = C(dfr[[colnm]],contrasts[[colnm]]);
            }
          } else
          if (!(is.factor(dfr[[colnm]])||is.numeric(dfr[[colnm]]))){
            #--convert column to factor if not already a factor or numeric
            levs = unique(as.character(dfr[[colnm]]));
            dfr[[colnm]] = factor(as.character(dfr[[colnm]]),levels=levs);
            if (is.list(contrasts)){
              if (colnm %in% names(contrasts)){
                if (!is.logical(contrasts[[colnm]])) dfr[[colnm]] = C(dfr[[colnm]],contrasts[[colnm]]);
              }
            }
          }
        }
      }
    }
    return(dfr);
}
# dfrp = convertToDimsFactors(dfrSparse); str(dfrp);
# dfrp = convertToDimsFactors(dfrp);      str(dfrp);
#'

#' @title Drop selected dimensions from a dimensions map
#' @description Function to drop selected dimensions from a dimensions map.
#' @param dms - dimensions map
#' @param drop - vector (character or integer) with dimensions to drop
#' @return dimensions map without dropped dimensions
#' @details If \code{drop} is an integer vector, the values
#' specify the order of dimensions to drop.
#'
#' @export
#'
dropDims<-function(dms,drop){
    dmnms = attr(dms,"dmnms");
    if (is.numeric(drop)){
      drop = dmnms[drop];
    }
    keep = dmnms[!(dmnms %in% drop)];
    return(keepDims(dms,keep));
}

#'
#' @title Keep only selected dimensions from a dimensions map
#' @description Function to keep only selected dimensions from a dimensions map.
#' @param dms - dimensions map
#' @param keep - vector (character or integer) with dimensions to keep
#' @return dimensions map with only kept dimensions
#' @details If \code{keep} is an integer vector, the values
#' specify the order of dimensions to keep.
#'
#' @import dplyr
#'
#' @export
#'
keepDims<-function(dms,keep){
    dmnms = attr(dms,"dmnms");
    dmlvs = attr(dms,"dmlvs");
    dmlns = attr(dms,"dmlns");
    if (is.numeric(keep)){
      keep = dmnms[keep];
    }
    c1 = names(dms)[1];
    kept = dms |> dplyr::select(dplyr::all_of(c(c1,keep))) |>
                  dplyr::distinct(dplyr::pick(dplyr::all_of(keep))) |>
                  dplyr::mutate(i=dplyr::row_number(),.before=1);
    names(kept)[1] = c1;
    attr(kept,"dmnms")<-keep;
    attr(kept,"dmlvs")<-dmlvs[keep];
    attr(kept,"dmlns")<-dmlns[keep];
    return(kept);
}

#'
#' @title Create a dims aggregator map
#' @description Function to create an aggregator map from one set of dimensions to another.
#'
#' @param dms_frm - aggregating from dimensions map
#' @param dms_to - aggregating to dimensions map
#' @param keepOrigDims - flag (T/F) to keep original dimensions (T) or aggregated dimensions (F)
#'
#' @return map from original dimensions to aggregated dimensions.
#'    * column "idx_to" gives index into vector with aggregated dimensions
#'    * column "idx_from" gives index into vector with original dimensions
#'
#' @import  dplyr
#'
#' @export
#'
createAggregatorMap<-function(dms_frm,dms_to,keepOrigDims=FALSE){
    nms_fr = attr(dms_frm,"dmnms");
    nms_to = attr(dms_to, "dmnms");
    names(dms_frm)[1] = "idx_from";
    names(dms_to)[1]  = "idx_to";
    agg = dms_frm |> dplyr::inner_join(dms_to,by=nms_to);
    if (keepOrigDims){
      agg = agg |> dplyr::relocate(idx_to,idx_from,.before=1);
    } else {
      agg = agg |> dplyr::select(idx_to,idx_from,dplyr::all_of(nms_to));
    }
    return(agg);
}
#aggMap = createAggregatorMap(dfrSparse,kept)

getSubset.DimsMap<-function(dm,...){
  cat("length =",...length(),"\n")
  if (...length()==1){
    y = ...elt(1);
    if (is.DimsMap(y)) {
      return(dm |> dplyr::inner_join(y |> dplyr::select(!1)));
    } else
    if (...names()[1]=="row") {return(dm[y,]);} else
    if (...names()[1]=="idx") {
      idx = names(x)[1];
      return(dm |> dplyr::filter(.data[[idx]]==y));
    } else {
      cat("evaluating single element dots list as DimsMap\n");
      y = createSparseDimsMap(rlang::list2(...));
      return(dm |> dplyr::inner_join(y |> dplyr::select(!1)));
    }
  } else {
    cat("evaluating multiple element dots list\n");
    y = createSparseDimsMap(...);
    return(dm |> dplyr::inner_join(y |> dplyr::select(!1)));
  }
  return(NULL);
}
