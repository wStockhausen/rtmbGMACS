#--extract natural mortality parameters
#'
#' @title Extract natural mortality parameters from parameter info list
#' @description Function to extract natural mortality parameters from parameter info list.
#' @param lst - parameter info list from [readParamInfo_NaturalMortality()]
#' @param dims - `dims` list from `setupModelDims`
#' @return a list (see details)
#' @details The format of the returned list depends on the `option` specified in `lst`.
#' If `lst$option` is "data" (i.e., the input parameter values are fixed values of weight-at-size
#' specified at least for every population category), the output list has elements
#' \itemize{
#' \item{option - format option}
#' \item{params - vector of weight-at-size values}
#' \item{map - vector for the RTMB `map` list indicating the parameters are fixed; i.e. not to be estimated}
#' \item{dfrIdx2Pars - dataframe with columns "pidx", dimension names (y,s,r,x,m,p,z), and parameter specs IV,LB,UB,phz,PriorType,Pr1,Pr2
#' indicating the un-expanded mapping from `pidx` to the parameter input values ("IV"). `pidx` is the index into the parameter vector.
#' The dimension columns are incidental and reflect the specification from the input text file.}
#' \item{dfrDims2Pars - dataframe with columns "pidx" and the dimension names (y,s,r,x,m,p,z).
#' The dimension levels are fully expanded such that each row corresponds to a single "multidimensional" level.
#' The value of `pidx` for each row indicates the index of the associated weight-at-size value in the parameter vector.}
#' }
#'
#' If `lst$option` is "function", the output list has elements
#' \itemize{
#' }
#'
#' @import dplyr
#'
#' @export
#'
extractParamInfo_NaturalMortality<-function(lst,
                                     dims=NULL,
                                     verbose=TRUE){
  if (verbose) message("starting extractParameters_NaturalMortality.")
  if (FALSE){
    #--NOTE: run this section if you are just stepping through the code for development purposes
    ##--assumes you've created `dims` via `dims = setupModelDims();`
    verbose = TRUE;
    lst = res;#--assumed output from `res  = readParamInfo_NaturalMortality(conn,verbose=FALSE);`
  }
  # #--create list of all dimension levels in `dims$dmsYSN` to convert "all"'s to pop levels----
  # ##--listAlls is a list with all individual dimension levels, by individual dimension y, s, r, x, m, a, p, z
  # lstAlls = NULL;
  # if (!is.null(dims)) lstAlls = alls_GetLevels(dims$dmsYSC,verbose=verbose);

  #--expand natural mortality parameter information----
  if (tolower(lst$option)=="function"){
    ##--option == "function"----
    ##--inputs are functions and parameters definitions
    out = extractParamInfoFunctionType1(lst,dims$dmsYSC,"NaturalMortality",verbose);

    # ####--functions----
    # if (verbose) message("in extractParameters_NaturalMortality: processing functions.")
    # dfrFcns = lst$Fcns$dfr |> expandDataframe(lstAlls,verbose=verbose) |>
    #             dplyr::select(y,s,r,x,m,p,z,fcn,fcn_idx);
    # dfrFcns_RefLvls = NULL;
    # if (!is.null(lst$Fcns$reflvls$dfr))
    #   dfrFcns_RefLvls = lst$Fcns$reflvls$dfr |> expandDataframe(lstAlls,verbose=verbose);
    #
    # ###--main parameters----
    # if (verbose) message("in extractParameters_NaturalMortality: processing main parameters.")
    # ####--mp_idx: input "index" for main parameters (not necessarily sequential)--use as join index for other tables
    # ####--add mpr_idx: sequential row "index" for main parameters parameter vector
    # ####--add in_mp_idx: sequential index "within" parameter names (??)
    # dfrMP1s = lst$MPs$dfr |>
    #             dplyr::mutate(mpr_idx=dplyr::row_number(),.after=mp_idx);      #--add index for vector of parameter values
    #             # dplyr::group_by(param) |>
    #             # dplyr::mutate(in_mp_idx=row_number(),.after=mpr_idx) |> #--add index within parameter name
    #             # dplyr::ungroup();
    # ####--transform initial and associated values
    # dfrMP2s = dfrMP1s  |> dplyr::rowwise() |>
    #             dplyr::mutate(IV=tf_apply(tform,IV),
    #                           LB=tf_apply(tform,LB),
    #                           UB=-tf_apply(tform,UB),
    #                           Pr1=tf_apply(tform,Pr1)) |>
    #             dplyr::ungroup();
    # # ####--reorder by index "within" parameter name
    # # dfrMP3s = dfrMP2s |> dplyr::arrange(param,pv_idx);
    # ####--assign initial values to parameters vector (indexed by `pv_idx`)
    # mpParams = dfrMP2s$IV;
    # ####--expand by all dimensions and drop values info
    # dfrMPs = dfrMP2s |> expandDataframe(lstAlls,verbose=verbose) |>
    #            dplyr::select(y,s,r,x,m,p,z,fcn_idx,grp_idx,param,mp_idx,mpr_idx);
    #
    # ###--offset parameters----
    # if (verbose) message("in extractParameters_NaturalMortality: processing offset parameters.")
    # dfrOPs = NULL;
    # if (!is.null(lst$OPs$dfr)){
    # ####--mp_idx: "index" identifying main parameters (not necessarily sequential)
    # ####--op_idx: "index" identifying offset parameters (not necessarily sequential)
    # ####--add opr_idx: sequential row "index" for offset parameters parameter vector
    # ####--add in_op_idx: sequential index "within" parameter names (??)
    #   dfrOP1s = lst$OPs$dfr |>
    #               dplyr::mutate(opr_idx=dplyr::row_number(),.after=op_idx);        #--add index for vector of parameter values
    #               # dplyr::group_by(param) |>
    #               # dplyr::mutate(in_off_idx=row_number(),.after=opr_idx) |> #--add offset index "within" parameter name
    #               # dplyr::ungroup() |>
    #   ####--transform initial and associated values
    #   dfrOP2s = dfrOP1s  |> dplyr::rowwise() |>
    #               dplyr::mutate(IV=tf_apply(tform,IV),
    #                             LB=tf_apply(tform,LB),
    #                             UB=-tf_apply(tform,UB),
    #                             Pr1=tf_apply(tform,Pr1))|>
    #               dplyr::ungroup();
    #   ####--reorder by offset index "within" parameter name
    #   # dfrOP3s = dfrOP2s |> dplyr::arrange(param,pv_idx);
    #   ####--assign initial values to parameters vector (indexed by `pv_idx`)
    #   opParams = dfrOP2s$IV;
    #   ####--expand by all dimensions and drop values info
    #   dfrOPs  = dfrOP2s |> expandDataframe(lstAlls,verbose=verbose) |>
    #               dplyr::select(y,s,r,x,m,p,z,param,mp_idx,op_idx,opr_idx,op_type);
    # }
    #
    # ###--devs parameter vectors----
    # if (verbose) message("in extractParameters_NaturalMortality: processing devs parameter vectors.")
    # dfrDPs = NULL; dfrDPs_RefLvls = NULL;
    # if (!is.null(lst$DPs$dfr)){
    #   ####--add devs vector index "within" parameter names
    #   dfrDP1s = lst$DPs$dfr; # |>
    #   ####--transform initial and associated values
    #   dfrDP2s = dfrDP1s  |> dplyr::rowwise() |>
    #               dplyr::mutate(IV=tf_apply(tform,IV),
    #                             LB=tf_apply(tform,LB),
    #                             UB=-tf_apply(tform,UB),
    #                             Pr1=tf_apply(tform,Pr1)) |>
    #               dplyr::ungroup();
    #   ####--create parameter index within each dev vector
    #   lstDP4s = list();
    #   for (rw in 1:nrow(dfrDP2s)){
    #     dfrTmp = dfrDP2s[rw,];
    #     xpnds    = stringr::str_split_1(dfrTmp$expand_across[1],",");
    #     lstTmp  = list();
    #     dfrTmpp = dfrTmp |> dplyr::select(!any_of(c("expand_across",xpnds)))
    #     for (xpnd in xpnds){
    #       sm_xpnd = rlang::sym(xpnd);
    #       if (dfrTmp[[xpnd]][1]=="all"){
    #         vals = as.character((lstAlls[[xpnd]] |> dplyr::filter({{sm_xpnd}}=="all"))$value);
    #       } else {
    #         vals = as.character(eval(parse(text=dfrTmp[[xpnd]][1])));
    #       }
    #       lstTmp[[xpnd]] = dfrTmpp  |>
    #                           dplyr::cross_join(tibble::tibble({{sm_xpnd}}:=vals));
    #       rm(sm_xpnd,vals);
    #     }
    #     lstDP4s[[rw]] = dplyr::bind_rows(lstTmp) |>
    #                       dplyr::mutate(in_dv_idx=dplyr::row_number(),.after=dv_idx);
    #     rm(dfrTmp,dfrTmpp,lstTmp,xpnd,xpnds);
    #   }##--rw loop
    #   dfrDP4s = dplyr::bind_rows(lstDP4s); rm(lstDP4s);
    #   ####--add dev parameter index
    #   dfrDP4s = dfrDP4s |> dplyr::mutate(dpr_idx=dplyr::row_number(),.after=in_dv_idx);
    #   ####--assign initial values to parameters across vectors (dfrDP4s$dev_par_idx identifies index into vector)
    #   dpParams = dfrDP4s$IV;
    #   dfrDPs  = dfrDP4s |> expandDataframe(lstAlls,verbose=verbose) |>
    #               dplyr::select(y,s,r,x,m,p,z,param,mp_idx,dv_idx,in_dv_idx,dpr_idx,dv_type);
    #   ####--get reference levels information
    #   if (!is.null(lst$DPs$reflvls$dfr))
    #     dfrDPs_RefLvls = lst$DPs$reflvls$dfr |> expandDataframe(lstAlls,verbose=verbose);
    # }
    #
    # ###--RE parameter vectors----
    # if (verbose) message("in extractParameters_NaturalMortality: processing RE parameter vectors.")
    # dfrREs = NULL; dfrDPs_RefLvls = NULL;
    # if (!is.null(lst$REs$dfr)){
    #   ####--add RE vector index "within" parameter names
    #   dfrRE1s = lst$REs$dfr;
    #   ####--transform initial and associated values
    #   dfrRE2s = dfrRE1s  |> dplyr::rowwise() |>
    #               dplyr::mutate(IV=tf_apply(tform,IV),
    #                             LB=tf_apply(tform,LB),
    #                             UB=-tf_apply(tform,UB),
    #                             Pr1=tf_apply(tform,Pr1)) |>
    #               dplyr::ungroup();
    #   ####--create parameter index within each RE vector
    #   lstRE4s = list();
    #   for (rw in 1:nrow(dfrRE2s)){
    #     dfrTmp = dfrRE2s[rw,];
    #     xpnds    = stringr::str_split_1(dfrTmp$expand_across[1],",");
    #     lstTmp  = list();
    #     dfrTmpp = dfrTmp |> dplyr::select(!any_of(c("expand_across",xpnds)))
    #     for (xpnd in xpnds){
    #       sm_xpnd = rlang::sym(xpnd);
    #       if (dfrTmp[[xpnd]][1]=="all"){
    #         vals = as.character((lstAlls[[xpnd]] |> dplyr::filter({{sm_xpnd}}=="all"))$value);
    #       } else {
    #         vals = as.character(eval(parse(text=dfrTmp[[xpnd]][1])));
    #       }
    #       lstTmp[[xpnd]] = dfrTmpp  |>
    #                           dplyr::cross_join(tibble::tibble({{sm_xpnd}}:=vals));
    #       rm(sm_xpnd,vals);
    #     }
    #     lstRE4s[[rw]] = dplyr::bind_rows(lstTmp) |>
    #                       dplyr::mutate(in_dv_idx=dplyr::row_number(),.after=dv_idx);
    #     rm(dfrTmp,dfrTmpp,lstTmp,xpnd,xpnds);
    #   }##--rw loop
    #   dfrRE4s = dplyr::bind_rows(lstRE4s); rm(lstRE4s);
    #   ####--add RE parameter index
    #   dfrRE4s = dfrRE4s |> dplyr::mutate(rpr_idx=dplyr::row_number(),.after=in_rv_idx);
    #   ####--assign initial values to parameters across vectors
    #   reParams = dfrRE4s$IV;
    #   dfrREs  = dfrRE4s |> expandDataframe(lstAlls,verbose=verbose) |>
    #               dplyr::select(y,s,r,x,m,p,z,param,mp_idx,rv_idx,in_rv_idx,rpr_idx,rv_type);
    #   ####--get reference levels information
    #   if (!is.null(lst$REs$reflvls$dfr))
    #     dfrREs_RefLvls = lst$REs$reflvls$dfr |> expandDataframe(lstAlls,verbose=verbose);
    # }
    #
    # ###--parameter-related environmental covariates----
    # if (verbose) message("in extractParameters_NaturalMortality: processing parameter-related environmental covariates.")
    # dfrPECs = NULL; dfrPECs_RefLvls = NULL;
    # if (!is.null(lst$PECs$dfr)){
    #   dfrPECs = lst$PECs$dfr |> expandDataframe(lstAlls,verbose=verbose);
    #   if (!is.null(lst$PECs$reflvls$dfr))
    #     dfrPECs_RefLvls = lst$PECs$reflvls$dfr |> expandDataframe(lstAlls,verbose=verbose);
    # }
    #
    # ###--functions-related environmental covariates----
    # if (verbose) message("in extractParameters_NaturalMortality: processing function-related environmental covariates.")
    # dfrFECs = NULL; dfrFECs_RefLvls = NULL;
    # if (!is.null(lst$FECs$dfr)){
    #   dfrFECs = lst$FECs$dfr |> expandDataframe(lstAlls,verbose=verbose);
    #   if (!is.null(lst$FECs$reflvls$dfr))
    #     dfrFECs_RefLvls = lst$FECs$reflvls$dfr |> expandDataframe(lstAlls,verbose=verbose);
    # }
    #
    # ###--functional priors----
    # if (verbose) message("in extractParameters_NaturalMortality: processing functional priors.")
    # dfrFPs = NULL; dfrFPs_RefLvls = NULL;
    # if (!is.null(lst$FPs$dfr)){
    #   dfrFPs = lst$FPs$dfr |> expandDataframe(lstAlls,verbose=verbose);
    # }
    #
    # ###--create full indexing dataframe----
    # ####--TODO: need to add:
    # ####--RE covariance info
    # ####--parameter-related covariates
    # ####--function-related covariates
    # dfrCmbs = dfrFcns |>
    #             dplyr::left_join(dfrMPs,by = dplyr::join_by(y, s, r, x, m, p, z, fcn_idx)) |>
    #             dplyr::mutate(idx=paste0(param," + ",mpr_idx));
    # if (!is.null(dfrOPs))
    #   dfrCmbs = dfrCmbs |>
    #               dplyr::left_join(dfrOPs,by = dplyr::join_by(y, s, r, x, m, p, z, mp_idx, param)) |>
    #               dplyr::mutate(idx=paste0(idx," + ",opr_idx));
    # if (!is.null(dfrDPs))
    #   dfrCmbs = dfrCmbs |>
    #               dplyr::left_join(dfrDPs,by = dplyr::join_by(y, s, r, x, m, p, z, mp_idx, param)) |>
    #               dplyr::mutate(idx=paste0(idx," + ",dpr_idx));
    # if (!is.null(dfrREs))
    #   dfrCmbs = dfrCmbs |>
    #               dplyr::left_join(dfrREs,by = dplyr::join_by(y, s, r, x, m, p, z, mp_idx, param)) |>
    #               dplyr::mutate(idx=paste0(idx," + ",rpr_idx));
    #
    # dfrUniqCmbs = dfrCmbs |>
    #                 dplyr::select(!c(y, s, r, x, m, p, z)) |>
    #                 dplyr::distinct();
    #
    # dfrHCs = dfrCmbs  |>
    #            dplyr::select(y,s,r,x,m,p,z,fcn,fcn_idx,grp_idx,param,idx) |>
    #            tidyr::pivot_wider(names_from=param,values_from=idx);
    #
    # dfrUHCs = dfrHCs |>
    #             dplyr::select(!c(y,s,r,x,m,p,z)) |>
    #             dplyr::distinct();
    #
    # ###--create output list----
    # dfrRefLvlsFcns=NULL;
    # if (!is.null(lst$Fcns$reflvls)) dfrRefLvlsFcns=lst$Fcns$reflvls$dfr;
    # dfrRefLvlsDPs=NULL;
    # if (!is.null(lst$DPs$reflvls)) dfrRefLvlsDPs=lst$DPs$reflvls$dfr;
    # dfrRefLvlsPECs=NULL;
    # if (!is.null(lst$PECs$reflvls)) dfrRefLvlsPECs=lst$PECs$reflvls$dfr;
    # dfrRefLvlsFECs=NULL;
    # if (!is.null(lst$FECs$reflvls)) dfrRefLvlsFECs=lst$FECs$reflvls$dfr;
    # out = list(option="function",
    #            Fcns=list(dfrIdxs=lst$Fcns$dfr,
    #                      dfrRefLvls=dfrRefLvlsFcns,
    #                      dfrDims2Idxs=dfrFcns),
    #            MPs=list(dfrMP1s=get0("dfrMP1s",ifnotfound=NULL),
    #                     dfrMP2s=get0("dfrMP2s",ifnotfound=NULL),
    #                     # dfrMP3s=get0("dfrMP3s",ifnotfound=NULL),
    #                     dfrMPs =get0("dfrMPs", ifnotfound=NULL),
    #                     params=get0("mpParams",ifnotfound=NULL)),
    #            OPs=list(dfrOP1s=get0("dfrOP1s",ifnotfound=NULL),
    #                     dfrOP2s=get0("dfrOP2s",ifnotfound=NULL),
    #                     # dfrOP3s=get0("dfrOP3s",ifnotfound=NULL),
    #                     dfrOPs =get0("dfrOPs", ifnotfound=NULL),
    #                     params=get0("opParams",ifnotfound=NULL)),
    #            DPs=list(dfrDP1s=get0("dfrDP1s",ifnotfound=NULL),
    #                     dfrDP2s=get0("dfrDP2s",ifnotfound=NULL),
    #                     # dfrDP3s=get0("dfrDP3s",ifnotfound=NULL),
    #                     dfrDP4s=get0("dfrDP4s",ifnotfound=NULL),
    #                     dfrDPs =get0("dfrDPs", ifnotfound=NULL),
    #                     params=get0("dpParams",ifnotfound=NULL)),
    #            REs=list(dfrRE1s=get0("dfrRE1s",ifnotfound=NULL),
    #                     dfrRE2s=get0("dfrRE2s",ifnotfound=NULL),
    #                     dfrRE4s=get0("dfrRE4s",ifnotfound=NULL),
    #                     dfrREs =get0("dfrREs", ifnotfound=NULL),
    #                     params=get0("reParams",ifnotfound=NULL)),
    #            PECs=list(dfrIdxs=lst$PECs$dfr,
    #                      dfrRefLvls=dfrRefLvlsPECs,
    #                      dfrDims2Idxs=dfrPECs,
    #                      params=lst$PECs$dfr$IV),
    #            FECs=list(dfrIdxs=lst$FECs$dfr,
    #                      dfrRefLvls=dfrRefLvlsFECs,
    #                      dfrDims2Idxs=dfrFECs,
    #                      params=lst$FECs$dfr$IV),
    #            dfrFPs=dfrFPs,
    #            dfrCmbs=dfrCmbs,
    #            dfrUniqCmbs=dfrUniqCmbs,
    #            dfrHCs=dfrHCs,
    #            dfrUHCs=dfrUHCs);
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
        message("in extractParamInfo_NaturalMortality: lst$dfr:")
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
    ###--add missing dimensions as "all"s
    dmnms=attr(dims$dmsYSC,"dmnms");
    for (dmnm in dmnms){
      if (!(dmnm %in% names(dfr))) dfr[[dmnm]] = "all";
    }
    dfr = dfr |> dplyr::select(pidx,y,s,r,x,m,p,z,IV,LB,UB,phz,PriorType,Pr1,Pr2);
    ###--extract parameter values----
    if (verbose) {
      message("in extractParameters_NaturalMortality: resolved 'data' option values");
      print(dfr);
    }
    pWatZ = dfr$IV;#--weight-at-size in kg
    map = list(pWatZ=factor(NA+pWatZ));#--NAs indicate fixed values
    if (verbose) message("in extractParameters_NaturalMortality: expanding dataframe.")
    ###--create list of all dimension levels in `dims$dmsYSN` to convert "all"'s to pop levels----
    ####--listAlls is a list with all individual dimension levels, by individual dimension y, s, r, x, m, a, p, z
    lstAlls = NULL;
    if (!is.null(dims)) lstAlls = alls_GetLevels(dims$dmsYSC,verbose=verbose);
    dfrp = dfr |>
              dplyr::select(!c(IV,LB,UB,phz,PriorType,Pr1,Pr2)) |>
              expandDataframe(lstAlls=lstAlls,verbose=verbose);#--expanded for "alls"
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
  #--test "data" option with "vertical" format-----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R/readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R/readParamInfo_NaturalMortality.R"))
  source(file.path(dirPrj,"R/extractParamInfo_NaturalMortality.R"))
  source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))
  dims = setupModelDims();
  conn = file.path(dirPrj,"testing/testNaturalMortality/inputSpecs_NaturalMortality.data-vertical.txt");
  res1v  = readParamInfo_NaturalMortality(conn,verbose=FALSE);
  res2v = extractParamInfo_NaturalMortality(res1v,dims,verbose=FALSE);
}

if (FALSE){
  #--test "data" option with "horizontal" format-----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R/readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R/readParamInfo_NaturalMortality.R"))
  source(file.path(dirPrj,"R/extractParamInfo_NaturalMortality.R"))
  source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))
  dims = setupModelDims();
  conn = file.path(dirPrj,"testing/testNaturalMortality/inputSpecs_NaturalMortality.data-horizontal.txt");
  res1h  = readParamInfo_NaturalMortality(conn,verbose=FALSE);
  res2h = extractParamInfo_NaturalMortality(res1h,dims,verbose=FALSE);
}

if (FALSE){
  #--test "function" option----
  dirPrj = rstudioapi::getActiveProject();
  source(file.path(dirPrj,"R/MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R/MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R/readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R/readParamInfo_NaturalMortality.R"))
  source(file.path(dirPrj,"R/extractParamInfo_NaturalMortality.R"))
  source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))
  dims = setupModelDims();
  conn = file.path(dirPrj,"testing/testNaturalMortality/inputSpecs_NaturalMortality.function.txt");
  res1f  = readParamInfo_NaturalMortality(conn,verbose=FALSE);
  res2f = extractParamInfo_NaturalMortality(res1f,dims,verbose=TRUE);
}
