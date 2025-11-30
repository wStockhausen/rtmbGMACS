readFile_CTL<-function(fn){
  #--read in CTL file----
  str = readr::read_lines(fn) |>
        paste("\n",collapse=""); #--need this as a single string

  #--process CTL file----
  strv = str |> splitText() |> extractLines("PROCESS_ALLOMETRY","END_ALLOMETRY");
  strv |> extractLines("CODE","END_CODE") |> evalTextAsCode(frame=.GlobalEnv);#--change to 0 in objfun

  ##--process function information----
  dfrFcns = strv |> extractLines("FUNCTIONS","END_FUNCTIONS")   |> evalTextAsDataframe();
  ###--expand function info to relevant model dimensions----
  lstFcnsExp = list();
  for (row in 1:nrow(dfrFcns)){
    rw = dfrFcns[row,];
    lstFcnsExp[[row]] = get(rw$frame) |>
                          dplyr::mutate(fcn_idx=rw$fcn_idx,fcn_id=rw$fcn_id,fcn_nm=rw$fcn_nm,
                                        params=list(stringr::str_split_1(rw$params,",")));
  }
  dfrFcnsExp = dplyr::bind_rows(lstFcnsExp);
  rm(row,rw,lstFcnsExp);

  ##--process parameters equation information
  dfrPEQs = strv |> extractLines("PARAM_EQS","END_PARAM_EQS") |> evalTextAsDataframe();
  par_ids = dfrPEQs$par_id;
  ###--fixed effects parameter values
  strFEVs = strv |> extractLines("FE_VALS","END_FE_VALS");
  lstFEVs = list();
  if (length(strFEVs)>0){
    for (par_id in par_ids){
      strp = strFEVs |> extractLines(par_id,paste0("END_",par_id));
      if (length(strp)>0)
        lstFEVs[[par_id]]  = strp |> evalTextAsDataframe() |> dplyr::mutate(par_id=par_id,.before=1) |> expand_mirrors();
      rm(strp);
    }
    rm(par_id);
  }#--FEVs
  ###--environmental covariates parameter values
  strECVs = strv |> extractLines("EC_VALS","END_EC_VALS");
  lstECVs = list();
  if (length(strECVs)>0){
    for (par_id in par_ids){
      strp = strECVs |> extractLines(par_id,paste0("END_",par_id));
      if (length(strp)>0)
        lstECVs[[par_id]]  = strp |> evalTextAsDataframe() |> dplyr::mutate(par_id=par_id,.before=1) |> expand_mirrors();
      rm(strp);
    }
    rm(par_id);
  }#--ECVs
  ###--random effects parameter values
  strREVs = strv |> extractLines("RE_VALS","END_RE_VALS");
  lstREVs = list();
  if (length(strREVs)>0){
    for (par_id in par_ids){
      strp = strREVs |> extractLines(par_id,paste0("END_",par_id));
      if (length(strp)>0)
        lstREVs[[par_id]]  = strp |> evalTextAsDataframe() |> dplyr::mutate(par_id=par_id,.before=1) |> expand_mirrors();
      rm(strp);
    }
    rm(par_id);
  }#--REVs
  dfrFcnParInfo = dfrFcns |> dplyr::full_join(dfrPEQs,by="fcn_idx");
  return(list(dfrFcnParInfo=dfrFcnParInfo,dfrFcnsExp=dfrFcnsExp,lstFEVs=lstFEVs,lstECVs=lstECVs,lstREVs=lstREVs))
}
