#--extract alllometry parameters
#'
#' @title Extract allometry parameters from parameter info list
#' @description Function to extract allometry parameters from parameter info list.
#' @param lst - parameter info list from [readParamInfo_Allometry()]
#' @param dims - `dims` list from `setupModelDims`
#' @return a list
#' @details TODO: fill in.
#' @export
#'
extractParamInfo_Allometry<-function(lst,
                                      dims=NULL,
                                      verbose=TRUE){
  if (verbose) message("starting extractParameters_Allometry.")
  if (FALSE){
    #--NOTE: run this section if you are just stepping through the code for development purposes
    ##--assumes you've created `dims` via `dims = setupModelDims();`
    verbose = TRUE;
    lst = res;#--assuemed output from `res  = readParamInfo_Allometry(conn,verbose=FALSE);`
  }
  #--create list of all dimension levels in `dims$dmsYSN` to convert "all"'s to pop levels----
  ##--listAlls is a list with all individual dimension levels, by individual dimension y, s, r, x, m, a, p, z
  lstAlls = NULL;
  if (!is.null(dims)) lstAlls = alls_GetLevels(dims$dmsYSN,verbose=verbose);

  #--create output list----
  out = list();

  #--expand allometry parameter information----
  if (tolower(lst$option)=="function"){
    ##--option == "function"----
    ##--inputs are functions and parameters definitions

    ####--functions----
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
    ####--combine parameter-related indices----
    dfrp = out$MPs$dfrDims2Idxs |>
             dplyr::select(!c(tform,IV,LB,UB,phz,PriorType,Pr1,Pr2)) |>
             dplyr::select(c(!dplyr::matches("par_idx|param"),param,par_idx)) |>
             dplyr::arrange(dplyr::pick(dplyr::any_of(c("y","r","x","m","a","p","z"))))
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
  } else if (tolower(lst$option)=="data"){
    ##--option == "data"----
    ###--inputs are data (fixed values)
    tform = stringr::str_starts(lst$transform,"tf_") |>
              ifelse(lst$transform,paste0("tf_",lst$transform));
    eval(parse(text=paste0("tf<-",tform)));#--evaluate transformation
    if (is.character(lst$dfr$value[1])){
      ###--values are character strings----
      ###--need to evaluate and transform `value`s to get (fixed) parameter values
      if (verbose){
        message("in extractParamInfo_Allometry: lst$dfr:")
        print(lst$dfr);
      }
      dfrp = lst$dfr |> dplyr::select(!c(z,value));
      dfrv = lst$dfr |> dplyr::select(vaue);
      lst1 = list();
      for (r in 1:nrow(lst$dfr)){
        lst1[[r]] =  dfrp[r,] |>
                   dplyr::cross_join(
                     tibble::as_tibble(eval(parse(text=dfrv[r,]))) |>
                       tidyr::pivot_longer(tidyselect::everything(),names_to="z",values_to="value")
                   )
      }
      dfr = dplyr::bind_rows(lst1);
      rm(dfrp,dfrv,lst1);
    } else {
      ####--values are numeric----
      #####--just copy lst4dfr
      dfr = lst$dfr;
    }
    ###--transform values, add bounds, phase < 0, priors info, parameter index----
    dfr = dfr |>
            dplyr::mutate(value=tf(value),
                          LB=-Inf,
                          UB= Inf,
                          phz=-1,
                          PriorType="none",
                          Pr1=0,
                          Pr2=999) |>
            dplyr::mutate(pidx=dplyr::row_number(),.before=1) |> #--parameter index
            dplyr::rename(IV=value);
    ###--extract parameter values----
    if (verbose) {
      message("in extractParameters_Allometry: resolved 'data' option values");
      print(dfr);
    }
    pWatZ = dfr$IV;#--weight-at-size in kg
    map = list(pWatZ=factor(NA+pWatZ));#--NAs indicate fixed values
    if (verbose) message("in extractParameters_Allometry: expanding dataframe.")
    dfrp = dfr |> expandDataframe(lstAlls=lstAlls,verbose=verbose);#--expanded for "alls"
    out = list(option="data",
               params=pWatZ,
               map=map,
               dfrIdx2Pars=dfr,
               dfrDims2Pars=dfrp);
  } else {
    ##--option not recognized----
    stop("option not recognized.")
  }
  return(out);
}

if (FALSE){
  #--test "data" option-----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R/readParamInfo_Allometry.R"))
  source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))
  conn = file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.data-vertical.txt");
  dims = setupModelDims();
  res  = readParamInfo_Allometry(conn,verbose=FALSE);
  res1 = extractParamInfo_Allometry(res,dims,verbose=FALSE);

}

if (FALSE){
  #--test "function" option----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R/readParamInfo_Allometry.R"))
  conn = file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.function.txt");
  dims = setupModelDims();
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
}
