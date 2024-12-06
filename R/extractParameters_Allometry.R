#--extract alllometry parameters
#'
#' @title Extract allometry parameters from parameter info list
#' @description Function to extract allometry parameters from parameter info list.
#' @param lst - parameter info list from [readParamInfo_Allometry()]
#' @param dims - dataframe with dimensions (columns) to expand
#' @return a list
#' @details TODO: fill in.
#' @export
#'
extractParameters_Allometry<-function(lst,
                                      dims=NULL,
                                      verbose=TRUE){
  if (verbose) message("starting extractParameters_Allometry.")
  if (FALSE){
    #TODO: delete this section--for development only when not running function
    dims = inputs$dims;
    lst = params_info$allometry;
  }
  lstAlls = NULL;
  if (!is.null(dims)) lstAlls = alls_GetLevels(dims,verbose=verbose);

  #--create output list----
  out = list();

  #--create `map` list for allometry----
  map = list();

  #--expand allometry parameter information
  if (lst$option=="function"){
    ##--inputs are functions and parameters definitions----
    ###--functions----
    if (verbose) message("in extractParameters_Allometry: processing functions.")
    dfrFcns = lst$Fcns$dfr |> expandDataframe(lstAlls,verbose=verbose);
    dfrFcns_RefLvls = NULL;
    if (!is.null(lst$Fcns$reflvls$dfr))
      dfrFcns_RefLvls = lst$Fcns$reflvls$dfr |> expandDataframe(lstAlls,verbose=verbose);

    ###--main parameters----
    if (verbose) message("in extractParameters_Allometry: processing main parameters.")
    dfrMPs = lst$MPs$dfr |> expandDataframe(lstAlls,verbose=verbose);

    ###--offset parameters----
    if (verbose) message("in extractParameters_Allometry: processing offset parameters.")
    dfrOPs = NULL;
    if (!is.null(lst$OPs$dfr))
      dfrOPs = lst$OPs$dfr |> expandDataframe(lstAlls,verbose=verbose);

    ###--devs parameters----
    if (verbose) message("in extractParameters_Allometry: processing devs parameters.")
    dfrDPs = NULL; dfrDPs_RefLvls = NULL;
    if (!is.null(lst$DPs$dfr)){
      dfrDPs = lst$DPs$dfr |> expandDataframe(lstAlls,verbose=verbose);
      if (!is.null(lst$DPs$reflvls$dfr))
        dfrDPs_RefLvls = lst$DPs$reflvls$dfr |> expandDataframe(lstAlls,verbose=verbose);
    }

    ###--parameter-related environmental covariates----
    if (verbose) message("in extractParameters_Allometry: processing parameter-related environmental covariates.")
    dfrPECs = NULL; dfrPECs_RefLvls = NULL;
    if (!is.null(lst$PECs$dfr)){
      dfrPECs = lst$PECs$dfr |> expandDataframe(lstAlls,verbose=verbose);
      if (!is.null(lst$PECs$reflvls$dfr))
        dfrPECs_RefLvls = lst$PECs$reflvls$dfr |> expandDataframe(lstAlls,verbose=verbose);
    }
    ###--functions-related environmental covariates----
    if (verbose) message("in extractParameters_Allometry: processing function-related environmental covariates.")
    dfrFECs = NULL; dfrFECs_RefLvls = NULL;
    if (!is.null(lst$FECs$dfr)){
      dfrFECs = lst$FECs$dfr |> expandDataframe(lstAlls,verbose=verbose);
      if (!is.null(lst$FECs$reflvls$dfr))
        dfrFECs_RefLvls = lst$FECs$reflvls$dfr |> expandDataframe(lstAlls,verbose=verbose);
    }

    ###--functional priors----
    if (verbose) message("in extractParameters_Allometry: processing functional priors.")
    dfrFPs = NULL; dfrFPs_RefLvls = NULL;
    if (!is.null(lst$FPs$dfr)){
      dfrFPs = lst$FPs$dfr |> expandDataframe(lstAlls,verbose=verbose);
    }
    ###--create output list----
    dfrRefLvlsFcns=NULL;
    if (!is.null(lst$Fcns$reflvls)) dfrRefLvlsFcns=lst$Fcns$reflvls$dfr;
    dfrRefLvlsDPs=NULL;
    if (!is.null(lst$DPs$reflvls)) dfrRefLvlsDPs=lst$DPs$reflvls$dfr;
    dfrRefLvlsPECs=NULL;
    if (!is.null(lst$PECs$reflvls)) dfrRefLvlsPECs=lst$PECs$reflvls$dfr;
    dfrRefLvlsFECs=NULL;
    if (!is.null(lst$FECs$reflvls)) dfrRefLvlsFECs=lst$FECs$reflvls$dfr;
    out = list(option="function",
               Fcns=list(dfrIdxs=lst$Fcns$dfr,
                         dfrRefLvls=dfrRefLvlsFcns,
                         dfrDims2Idxs=dfrFcns),
               MPs=list(dfrIdxs=lst$MPs$dfr,
                        dfrDims2Idxs=dfrMPs,
                        params=(lst$MPs$dfr |> dplyr::rowwise() |>
                                dplyr::mutate(IV=eval(parse(text=paste0(ifelse(stringr::str_starts(tform,"tf_"),
                                                                               tform,paste0("tf_",tform)),
                                                                        "(",IV,")")))))$IV),
               OPs=list(dfrIdxs=lst$OPs$dfr,
                        dfrDims2Idxs=dfrOPs,
                        params=(lst$OPs$dfr |> dplyr::rowwise() |>
                                dplyr::mutate(IV=eval(parse(text=paste0(ifelse(stringr::str_starts(tform,"tf_"),
                                                                               tform,paste0("tf_",tform)),
                                                                        "(",IV,")")))))$IV),
               DPs=list(dfrIdxs=lst$DPs$dfr,
                        dfrRefLvls=dfrRefLvlsDPs,
                        dfrDims2Idxs=dfrDPs,
                        params=(lst$DPs$dfr |> dplyr::rowwise() |>
                                dplyr::mutate(IV=eval(parse(text=paste0(ifelse(stringr::str_starts(tform,"tf_"),
                                                                               tform,paste0("tf_",tform)),
                                                                        "(",IV,")")))))$IV),
               PECs=list(dfrIdxs=lst$PECs$dfr,
                         dfrRefLvls=dfrRefLvlsPECs,
                         dfrDims2Idxs=dfrPECs,
                         params=lst$PECs$dfr$IV),
               FECs=list(dfrIdxs=lst$FECs$dfr,
                         dfrRefLvls=dfrRefLvlsFECs,
                         dfrDims2Idxs=dfrFECs,
                         params=lst$FECs$dfr$IV),
               dfrFPs=dfrFPs)
    ###--combine parameter-related indices----
    dfrp = out$MPs$dfrDims2Idxs |>
             dplyr::select(!c(tform,IV,LB,UB,phz,PriorType,Pr1,Pr2)) |>
             dplyr::select(c(!dplyr::matches("par_idx|param"),param,par_idx)) |>
             dplyr::arrange(dplyr::pick(dplyr::any_of(c("y","r","x","m","a","z"))))
    if (verbose) head(dfrp);
    if (!is.null(out$OPs$dfrDims2Idxs)){
      dfrp = dfrp |>
            dplyr::left_join(out$OPs$dfrDims2Idxs |> dplyr::select(!c(param,offset_type,tform,IV,LB,UB,phz,PriorType,Pr1,Pr2)));
      if (verbose) head(dfrp |> dplyr::filter(!is.na(off_idx)));
    }
    if (!is.null(out$DPs$dfrDims2Idxs)){
      dfrp = dfrp |>
               dplyr::left_join(out$DPs$dfrDims2Idxs |> dplyr::select(!c(param,dev_type,rw_type,`RE?`,expand_over,tform,IV,LB,UB,phz,PriorType,Pr1,Pr2)));
      if (verbose) head(dfrp |> dplyr::filter(!is.na(dev_idx)));
    }
    if (!is.null(out$PECs$dfrDims2Idxs)){
      dfrp = dfrp |>
                dplyr::nest_join(out$PECs$dfrDims2Idxs |> dplyr::select(!c(param,link_type,IV,LB,UB,phz,PriorType,Pr1,Pr2)),name="pcov_idxs");
      if (verbose) head(dfrp);
    }
    out$dfrDims2ParIdxs = dfrp;
  } else if (lst$option=="data"){
    ##--inputs are fixed values----
    tform = stringr::str_starts(lst$transform,"tf_") |>
              ifelse(lst$transform,paste0("tf_",lst$transform));
    eval(parse(text=paste0("tf<-",tform)));#--evaluate transformation
    if (is.character(lst$dfr$IV[1])){
      #--need to evaluate `IV`s to get (fixed) parameter values
      if (verbose){
        message("in extractParameters_Allometry: lst$dfr:")
        print(lst$dfr);
      }
      dfrp = lst$dfr |> dplyr::select(!c(z,IV));
      dfrv = lst$dfr |> dplyr::select(IV);
      lst1 = list();
      for (r in 1:nrow(lst$dfr)){
        lst1[[r]] =  dfrp[r,] |>
                   dplyr::cross_join(
                     tibble::as_tibble(eval(parse(text=dfrv[r,]))) |>
                       tidyr::pivot_longer(tidyselect::everything(),names_to="z",values_to="IV")
                   )
      }
      dfr = dplyr::bind_rows(lst1);
      rm(dfrp,dfrv,lst1);
    } else {
      dfr = lst$dfr;
    }
    ###--transform values, add bounds, phase < 0, priors info, parameter index----
    dfr = dfr |>
            dplyr::mutate(IV=tf(IV),
                          LB=-Inf,
                          UB= Inf,
                          phz=-1,
                          PriorType="none",
                          Pr1=0,
                          Pr2=999) |>
            dplyr::mutate(pidx=dplyr::row_number(),.before=1); #--parameter index
    ###--extract parameter values----
    if (verbose) {
      message("in extractParameters_Allometry: resolved 'data' option values");
      print(dfr);
    }
    pWatZ = dfr$IV;#--weight-at-size in kg
    map[["pWatZ"]] = factor(NA+pWatZ);#--NAs indicate fixed values
    if (verbose) message("in extractParameters_Allometry: expanding dataframe.")
    dfrp = dfr |> expandDataframe(lstAlls=lstAlls,verbose=verbose);#--expanded for "alls"
    out = list(option="data",
               params=pWatZ,
               map=map,
               dfrIdx2Pars=dfr,
               dfrDims2Pars=dfrp);
  } else {
    stop("option not recognized.")
  }
  return(out);
}

# #--"data" option----------------
# conn = file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.data.txt");
# res  = readParamInfo_Allometry(conn,verbose=FALSE);
# res1 = extractParameters_Allometry(res,dims$dms_yrxmaz,verbose=FALSE);
#
#--"function" option-----------
conn = file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.function.txt");
res  = readParamInfo_Allometry(conn,verbose=FALSE);
res1 = extractParameters_Allometry(res,dims$dms_yrxmaz,verbose=FALSE);

# dfrp = res1$MPs$dfrDims2Idxs |>
#          dplyr::select(!c(tform,IV,LB,UB,phz,PriorType,Pr1,Pr2)) |>
#          dplyr::select(c(!dplyr::matches("par_idx|param"),param,par_idx)) |>
#          dplyr::arrange(dplyr::pick(dplyr::any_of(c("y","r","x","m","a","z"))))
# head(dfrp);
# if (!is.null(res1$OPs$dfrDims2Idxs)){
#   dfrp = dfrp |>
#         dplyr::left_join(res1$OPs$dfrDims2Idxs |> dplyr::select(!c(param,offset_type,tform,IV,LB,UB,phz,PriorType,Pr1,Pr2)));
#   head(dfrp |> dplyr::filter(!is.na(off_idx)))
# }
# if (!is.null(res1$DPs$dfrDims2Idxs)){
#   dfrp = dfrp |>
#            dplyr::left_join(res1$DPs$dfrDims2Idxs |> dplyr::select(!c(param,dev_type,rw_type,`RE?`,expand_over,tform,IV,LB,UB,phz,PriorType,Pr1,Pr2)));
#   head(dfrp |> dplyr::filter(!is.na(dev_idx)))
# }
# if (!is.null(res1$PECs$dfrDims2Idxs)){
#   dfrp1 = dfrp |>
#             dplyr::nest_join(res1$PECs$dfrDims2Idxs |> dplyr::select(!c(param,link_type,IV,LB,UB,phz,PriorType,Pr1,Pr2)),name="pcov_idxs");
#
# }
#
# dfrf = res1$dfrFcns |> dplyr::select(!fcn) |> dplyr::select(c(!tidyselect::matches("fcn_idx"),fcn_idx));
# head(dfrf);
# if (!is.null(res1$dfrFECs)){
#   dfrf = dfrf |>
#         dplyr::nest_join(res1$dfrFECs |> dplyr::select(!c(fcn,link_type,IV,LB,UB,phz,PriorType,Pr1,Pr2)),name="fcov_idxs");
# }
# head(dfrf);
