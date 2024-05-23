#'
#' @title Convert dimension values to factor levels in a dimensions map
#' @description Function to convert dimension values to factor levels in a dimensions map.
#' @param dfr - dimensions map dataframe/tibble (e.g. output from [createSparseDimsMap()])
#'
#' @return - dimensions map with dimension values converted to factors
#'
#' @details Dimension columns that are already factors are not converted.
#'
#' @export
#'
createDimsFactors<-function(dfr){
    dmnms = attr(dfr,"dmnms");
    dmlvs = attr(dfr,"dmlvs");
    for (dmnm in dmnms){
        if (!is.factor(dfr[[dmnm]]))
            dfr[[dmnm]] = factor(dfr[[dmnm]],levels=dmlvs[[dmnm]]);
    }
    return(dfr);
}

# dfrp = createDimsFactors(dfrSparse); str(dfrp);
# dfrp = createDimsFactors(dfrp);      str(dfrp);

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
#' @param drop - vector (character or integer) with dimensions to keep
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

