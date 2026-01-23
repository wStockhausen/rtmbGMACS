readFile_CTLrev<-function(fn){
  #--read in CTL file----
  str = readr::read_lines(fn) |>
        paste("\n",collapse=""); #--need this as a single string

  #--process CTL file----
  strv = str |> splitText() |> extractLines("PROCESS_ALLOMETRY","END_ALLOMETRY");

  ##--process code and create objects----
  strv |> extractLines("CODE","END_CODE") |> evalTextAsCode(frame=.GlobalEnv);

  ##--process parameters equation information----
  dfrPEQs = strv |>
              extractLines("PARAM_EQS","END_PARAM_EQS") |>
              evalTextAsDataframe() |>
              dplyr::mutate(par_key=paste0(par_id,"_",par_idx),.before=1) |>
              dplyr::rename(par_frame=frame);
  par_keys = dfrPEQs$par_key;
  ###--fixed effects parameter values----
  strFEVs = strv |> extractLines("FE_VALS","END_FE_VALS");
  dfrFEVs = NULL;
  if (length(strFEVs) > 0) {
    dfrFEVs = strFEVs |>
                evalTextAsDataframe() |>
                dplyr::mutate(par_key=paste0(par_id,"_",par_idx)) |>
                expand_mirrors(); #--TODO: need to revise `expand_mirrors`
  }#--FEVs
  ###--environmental covariates parameter values----
  strECVs = strv |> extractLines("EC_VALS","END_EC_VALS");
  dfrECVs = NULL;
  if (length(strECVs)>0){
    dfrECVs = strECVs |>
                evalTextAsDataframe() |>
                dplyr::mutate(par_key=paste0(par_id,"_",par_idx)) |>
                expand_mirrors(); #--TODO: need to revise `expand_mirrors`
  }#--ECVs
  ###--random effects parameter values----
  strREVs = strv |> extractLines("RE_VALS","END_RE_VALS");
  dfrREVs = NULL;
  if (length(strREVs)>0){
    dfrREVs = strREVs |>
                evalTextAsDataframe() |>
                dplyr::mutate(par_key=paste0(par_id,"_",par_idx)) |>
                expand_mirrors(); #--TODO: need to revise `expand_mirrors`
  }#--REVs

  ##--process function information----
  dfrFcns = strv |> extractLines("FUNCTIONS","END_FUNCTIONS") |> evalTextAsDataframe()|>
              dplyr::rename(fcn_frame=frame) |>
              dplyr::rowwise() |>
              dplyr::mutate(par_keys = list(tibble::enframe(eval(parse(text=par_idxs))) |>
                                               dplyr::mutate(par_key=paste0(name,"_",value)) |>
                                               dplyr::select(par_key))) |>
              dplyr::ungroup();
  ###--map parameters to functions
  dfrMapParsToFcns = dfrFcns |> dplyr::select(fcn_idx,par_keys);
  lstMapParsToFcns = list();
  for (row in 1:nrow(dfrMapParsToFcns)){
    rw = dfrMapParsToFcns[row,];
    lstMapParsToFcns[[row]] = rw |> dplyr::select(fcn_idx) |>
                               dplyr::cross_join(rw$par_keys[[1]])
  }
  dfrMapParsToFcns = dplyr::bind_rows(lstMapParsToFcns);

  ###--expand function info to relevant model dimensions----
  lstFcnsExp = list();
  for (row in 1:nrow(dfrFcns)){
    #--testing: row = 2;
    rw = dfrFcns[row,];
    lstFcnsExp[[row]] = get(rw$fcn_frame) |>
                          dplyr::mutate(fcn_idx=rw$fcn_idx,fcn_id=rw$fcn_id,fcn_nm=rw$fcn_nm,
                                        vars=rw$vars,par_keys=rw$par_keys);
  }
  dfrFcnsExp = dplyr::bind_rows(lstFcnsExp) |>
                 dplyr::rowwise() |>
                 dplyr::mutate(par_keys=list(tibble::deframe(par_keys))) |> #--want to convert to character vector?
                 dplyr::ungroup();
  rm(row,rw,lstFcnsExp);

  return(list(dfrPEQs=dfrPEQs,dfrFEVs=dfrFEVs,dfrECVs=dfrECVs,dfrREVs=dfrREVs,
              dfrFcns=dfrFcns,dfrMapParsToFcns=dfrMapParsToFcns,dfrFcnsExp=dfrFcnsExp))
}
