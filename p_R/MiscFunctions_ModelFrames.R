#' @title Expand a functions dataframe to a model frame
#' @description Function to expand a functions dataframe to a model frame.
#' @param dfr - the functions dataframe
#' @return a `DimsMap` model frame for the functions
#'
#' @importFrom dplyr anti_join
#' @importFrom dplyr cross_join
#' @importFrom dplyr bind_rows
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @export
#'
expandFunctionsToModelFrame<-function(dfr){
  mfs = unique(dfr$frame);
  #--expand to first frame
  mfr = get(mfs[1]) |> dplyr::cross_join(dfr |> dplyr::filter(frame==mfs[1]));
  if (length(mfs)>1){
    for (i in 2:length(mfs)){
      #--testing: i=2;
      mfrp = get(mfs[i]) |> dplyr::cross_join(dfr |> dplyr::filter(frame==mfs[i]));
      mfr = dplyr::bind_rows(mfr |> dplyr::anti_join(mfrp |> dplyr::select(sparse_idx)),mfrp)
    }
  }
  return(mfr);
}

#--??
expandParamsToModelFrame<-function(dfr){
  bys = whichDims(names(dfr));
  mfs = unique(dfr$frame);
  #--expand to first frame
  mfr = get(mfs[1]) |> dplyr::left_join(dfr |> dplyr::filter(frame==mfs[1]),
                                        by=bys,relationship="many-to-many");
  if (length(mfs)>1){
    for (i in 2:length(mfs)){
      #--testing: i=2;
      mfrp = get(mfs[i]) |> dplyr::cross_join(dfr |> dplyr::filter(frame==mfs[i]));
      mfr = dplyr::bind_rows(mfr |> dplyr::anti_join(mfrp |> dplyr::select(sparse_idx)),mfrp)
    }
  }
  return(mfr);
}
