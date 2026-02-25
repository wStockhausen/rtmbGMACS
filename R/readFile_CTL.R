
#' @title Read CTL file TODO: finish this!
#' @description Function to read the CTL file.
#' @param fn - path/name to the CTL file
#' @return very long list
#'
#' @importFrom dplyr filter rename
#' @importFrom stringr fixed str_replace
#'
#' @export
#'
readFile_CTL<-function(fn){
  if (fn==""){
    dirPrj = rstudioapi::getActiveProject();
    dirThs = file.path(dirPrj,"testing/newTest_FilesCTL-Allometry");
    fn = file.path(dirThs,"CTL_Allometry_FormulaApproach.txt");
  }
  #--read in CTL file----
  str = readr::read_lines(fn) |>
        paste("\n",collapse=""); #--need this as a single string

  #--process CTL file----
  strv = str |> splitText();

  ##--process setup code and create objects----
  strv |> extractLines("SETUP_CODE","END_SETUP_CODE") |> evalTextAsCode(frame=.GlobalEnv);

  ##--parse allometry-specific code
  strv1 = strv |> extractLines("PROCESS_ALLOMETRY","END_ALLOMETRY");
  lstAllometryInfo = parseText_ProcessInfo(strv1);

  ##--parse code for other processes----

  return(namedList(lstAllometryInfo))
}
#--testing: lstCTL = readFile_CTL("");


parseText_ProcessInfo<-function(strv1){
  strv1 |> extractLines("CODE","END_CODE") |> evalTextAsCode(frame=.GlobalEnv);

  ###--process parameters equation information----
  dfrPEQs = strv1 |>
              extractLines("PARAM_EQS","END_PARAM_EQS") |>
              evalTextAsDataframe() |>
              # dplyr::mutate(par_key=paste0(par_id,"_",par_idx),.before=1) |>
              dplyr::rename(par_frame=frame);
  par_ids = dfrPEQs$par_id;
  call_strb = "createParamInfo(fm&,mf&,cs&,du&)";
  lstPIs = list();
  for (par_id_ in par_ids){
    #--testing: par_id_ = par_ids[1];
    rw = dfrPEQs |> dplyr::filter(par_id==par_id_);
    call_str = call_strb |> stringr::str_replace(stringr::fixed("fm&"),rw$formula) |>
                            stringr::str_replace(stringr::fixed("mf&"),rw$par_frame) |>
                            stringr::str_replace(stringr::fixed("cs&"),rw$contrasts) |>
                            stringr::str_replace(stringr::fixed("du&"),"TRUE");
    lstPIs[[par_id_]] = eval(str2expression(call_str));
  }

  ###--process parameter values information
  strv_pvs = strv1 |> extractLines("PARAM_VALS","END_PARAM_VALS");
  lstParVals = list();
  for (par_id_ in par_ids){
    #--testing: par_id_ = par_ids[1];
    strv_pvsp = strv_pvs |> extractLines(par_id_,paste0("END_",par_id_));
    strv_pvspp = strv_pvsp |> removeCommentLines();
    #--FEs
    ctr = 1;
    nFEs = extractNumericValue(strv_pvspp[ctr]); ctr<-ctr+1;
    dfrValInfo_FEs = NULL;
    if (nFEs>0) {
      dfrValInfo_FEs = extractDataframe(strv_pvspp,nFEs,ctr) |> dplyr::mutate(param=par_id_,.before=1);
      ctr = ctr+nFEs+1;
      }
    #--REs
    nReTerms = extractNumericValue(strv_pvspp[ctr]); ctr<-ctr+1;
    lstREs=NULL;
    if (nReTerms>0){
      piREs = lstPIs[[par_id_]]$piREs;
      lstREs = list();
      for (trm in piREs$reTerms){
        #--testing: trm = piREs$reTerms[1];
        #--TODO: need to read RE term here?
        phs = extractNumericValue(strv_pvspp[ctr]);  ctr<-ctr+1;
        nREs = extractNumericValue(strv_pvspp[ctr]); ctr<-ctr+1;
        dfrValInfo_REs = extractDataframe(strv_pvspp,nREs,ctr); ctr = ctr+nREs+1;
        cov_type = extractTextValue(strv_pvspp[ctr]);    ctr=ctr+1;
        nCovPars = extractNumericValue(strv_pvspp[ctr]); ctr = ctr+1;
        dfrValInfo_CPs = extractDataframe(strv_pvspp,nCovPars,ctr); ctr = ctr+nCovPars+1;
        lstREs[[trm]] = list(phs=phs,dfrValInfo=dfrValInfo_REs,lbls=piREs$re_pvs_lbls[[trm]],
                             cov_type=cov_type,nCovPars=nCovPars,dfrValInfo_CPs=dfrValInfo_CPs);
      }
    }
    #--SMs
    nSmTerms = extractNumericValue(strv_pvspp[ctr]); ctr<-ctr+1;
    lstSMs = NULL;
    if (nSmTerms>0){
      #--TODO: fill this in!
    }
    lstParVals[[par_id_]] = list(dfrValInfo_FEs=dfrValInfo_FEs,lstREs=lstREs,lstSMs=lstSMs);
  }#--par_ids loop

  lstParsInfo = list(dfrPEQs=dfrPEQs,lstPIs=lstPIs,lstParVals=lstParVals);

  ###--process function info
  dfrFEQs = strv1 |>
              extractLines("FUNCTIONS","END_FUNCTIONS") |>
              evalTextAsDataframe() |>
              dplyr::rename(fcn_frame=frame);
  fcn_ids = dfrFEQs$fcn_id;
  call_strb = "createFunctionInfo('fi&','fn&',mf&,'vr&',par_ids)";
  lstFIs = list();
  for (fcn_id_ in fcn_ids){
    #--testing: fcn_id_ = fcn_ids[1];
    rw = dfrFEQs |> dplyr::filter(fcn_id==fcn_id_);
    par_ids = stringr::str_split_1(rw$par_ids,",");
    call_str = call_strb |> stringr::str_replace(stringr::fixed("fi&"),rw$fcn_id) |>
                            stringr::str_replace(stringr::fixed("fn&"),rw$fcn) |>
                            stringr::str_replace(stringr::fixed("mf&"),rw$fcn_frame) |>
                            stringr::str_replace(stringr::fixed("vr&"),rw$vars);
    #cat(call_str,"\n");
    lstFIs[[fcn_id_]] = eval(str2expression(call_str));
  } #--fcn_ids loop
  lstFcnInfo = list(dfrFEQs=dfrFEQs,lstFIs=lstFIs);

  lstInfo = list(lstFcnInfo=lstFcnInfo,lstParsInfo=lstParsInfo);
  return(lstInfo);
}#--parseText_ProcessInfo
