extractParameters_InitialN<-function(lst,
                                     dims){
  if (FALSE){
    #TODO: delete this section--for development only when not running function
    dims=inputs$dims;
    lst = params_info$initN;
  }
  lstAlls = NULL;
  if (!is.null(dims)) lstAlls = alls_GetLevels(dims);

  #--expand initialN info----
  if (lst$option=="zero_pop") {
    ##--option=="zero_pop"
  } else if (lst$option=="steady-state_unfished") {
    ##--option=="steady-state_unfished"
  } else if (lst$option=="steady-state_fished") {
    ##--option=="steady-state_fished"
    pLnRini
  } else if (lst$option=="free_parameters_v1") {
    ##--option=="free_parameters_v1"
    dfrLnN0 = expandDataframe(lst$dfrLnN0) |>
                convertColsToNum(c("IV","LB","UB","phz","Pr1","Pr2"));
    pLnN0 = dfr$IV;
  } else if (lst$option=="free_parameters_v2") {
    ##--option=="free_parameters_v2"
    dfrLnTotInitPopSize = expandDataframe(lst$dfrLnTotInitPopSize) |>
                convertColsToNum(c("IV","LB","UB","phz","Pr1","Pr2"));
    pLnTotInitPopSize = dfrLnTotInitPopSize$IV;
    ref  = lst$dfrRef |> dplyr::mutate(z=as.character(z)) |>
                 dplyr::inner_join(dfrLnN0);
    #--TODO: set phase of refIdx to < 0 and
    #--      subtract IV[refIdx] from IV, LB, UB, and Pr1 for all rows
    dfrLnN0 = dfrLnN0 |> dplyr::mutate(IV=IV-ref$IV,
                                       LB=LB-ref$IV,
                                       UB=UB-ref$IV,
                                       Pr1=Pr1-ref$IV);
    dfrLnN0$phz[ref$idx] = -100;#--turn off estimation of reference
    pLnN0 = dfr$IV;
    return(list(option=lst$option,
                dfrLnTotInitPopSize=dfrLnTotInitPopSize,pLnTotInitPopSize=pLnTotInitPopSize,
                dfrLnN0=dfrLnN0,pLnN0=pLnN0));
  } else{
    stop(paste0("Initial N option '",lst$option,"' not recognized."));
  }
  return(NULL);
}




