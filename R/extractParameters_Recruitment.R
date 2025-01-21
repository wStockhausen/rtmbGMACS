extractParameters_Recruitment<-function(lst,
                                        dims=NULL){
  if (FALSE){
    #--TODO: for testing--for development only when not running function
    dims=inputs$dims
    lst = param_info$recruitment
  }
  #--get dims levels
  lstAlls = NULL;
  if (!is.null(dims)) lstAlls = alls_GetLevels(dims);

  #--expand parameter info time blocks
  dfrTBs = lst$lstTBs$dfrTBs |>
            dplyr::rowwise() |>
            dplyr::mutate(time_block=list(eval(parse(text=time_block)))) |>
            dplyr::ungroup() |>
            tidyr::unnest_longer(time_block,values_to="year") |>
            dplyr::select(time_block=name,year,label);
  #--expand process info
  dfrPrcs = expandDataframe(lst$lstPrcs$dfr,lstAlls[c("r","x","m","a","z")]);
  #--expand total recruitment specification
  dfrTR = expandDataframe(lst$lstTotRec,lstAlls[c("r","x","m","a","z")]);
  #--expand sex ratio specification
  #--expand size distribution specification

}


