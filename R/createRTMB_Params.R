#' @title Create the parameters information list for RTMB
#' @description Function to create the parameters information list for RTMB.
#' @param - list from [readFile_CTL()]
#' @return list
#' @details returned list has elements
#' \itemize{
#'  \item{lstRTMB_Params - list with initial parameter values}
#'  \item{lstRTMB_TForms - list with parameter value-level transform function names}
#'  \item{lstRTMB_LBs - list with parameter value lower bounds}
#'  \item{lstRTMB_UBs - list with parameter value upper bounds}
#'  \item{lstRTMB_Phases - list with parameter value estimation phases}
#'  \item{lstRTMB_Mirrors - list with parameter values mirroring}
#'  \item{lstRTMB_Xs - list with (sparse) design matrices}
#'  \item{lstRTMB_MMs - list with model matrices (as tibbles)}
#'  \item{lstRTMB_Map - RTMB "map" list}
#'  \item{dfrMMs - tibble with model matrix for all parameters expanded to the largest extent present}
#'  \item{vLinks - vector with parameter-level link function names}
#' }
#' @importFrom dplyr select
#' @importFrom tidyselect any_of
#' @export
#'
createRTMB_Params<-function(lst){
  #--extract parameter value-level information
  lstRTMB_Params  = list(); #--output list for RTMB parameters
  lstRTMB_TForms  = list(); #--output list for transforms from input parameter scales to RTMB scales
  lstRTMB_LBs     = list(); #--output list for RTMB parameter lower
  lstRTMB_UBs     = list(); #--output list for RTMB parameter upper bounds
  lstRTMB_Phases  = list(); #--output list for RTMB parameter phases
  lstRTMB_Mirrors = list(); #--output list for RTMB parameter mirrors
  lstRTMB_Xs      = list(); #--output list for RTMB design matrices
  lstRTMB_MMs     = list(); #--output list for RTMB model matrices
  lstRTMB_Map     = list(); #--output list for RTMB parameters "map" list
  lstMM = list();
  vREs = vector("character");
  par_ids = lst$dfrPEQs$par_id;
  for (par_id in par_ids){
    #--testing: par_id = par_ids[1];
    dfrMM = NULL;
    lstPIs = lst$lstPIs[[par_id]];
    lstPVs = lst$lstParVals[[par_id]];
    if (!is.null(lstPIs$piFEs)){
      #--fixed effects: want IVs, LBs, UBs, phase info, mirroring info, design matrix
      str = paste0("fe_",par_id);
      pars = tform(lstPVs$dfrValInfo_FEs$tform,lstPVs$dfrValInfo_FEs$IV);
      npars = length(pars);
      lstRTMB_Params[[str]]  = pars;
      lstRTMB_TForms[[str]]  = lstPVs$dfrValInfo_FEs$tform;
      lstRTMB_LBs[[str]]     = tform(lstPVs$dfrValInfo_FEs$tform,as.numeric(lstPVs$dfrValInfo_FEs$LB));
      lstRTMB_UBs[[str]]     = tform(lstPVs$dfrValInfo_FEs$tform,as.numeric(lstPVs$dfrValInfo_FEs$UB));
      lstRTMB_Phases[[str]]  = lstPVs$dfrValInfo_FEs$Phz;
      lstRTMB_Mirrors[[str]] = lstPVs$dfrValInfo_FEs$mirror;
      lstRTMB_Xs[[str]]      = lstPIs$piFEs$X;               #--design matrix (sparse matrix)
      lstRTMB_MMs[[str]]     = lstPIs$piFEs$model_matrix;    #--model matrix
      #--create "map" entry
      fac = 1:npars;
      #--get indices of mirrored parameters
      idxm = lstRTMB_Mirrors[[str]]>0;
      if (any(idxm)){
        #--set indices of mirrored parameters to those of the original parameters
        fac[idxm] = lstRTMB_Mirrors[[str]][idxm];
        #-- adjust lstRTMB_* entries for mirroring EXCEPT for lstRTMB_Xs, lstRTMB_MMs
        lstRTMB_Params[[str]]  = lstRTMB_Params[[str]][fac];
        lstRTMB_TForms[[str]]  = lstRTMB_TForms[[str]][fac];
        lstRTMB_LBs[[str]]     = lstRTMB_LBs[[str]][fac];
        lstRTMB_UBs[[str]]     = lstRTMB_UBs[[str]][fac];
        lstRTMB_Phases[[str]]  = lstRTMB_Phases[[str]][fac];
      }
      #--turn "off" estimation of parameter values with phase < 0
      idxp = lstRTMB_Phases[[str]]<0;#--indices of parameter values with phase < 0
      fac[idxp] = NA;                #--if "original"s are fixed, so are mirrored values
      rm(idxm,idxp);
      lstRTMB_Map[[str]] = factor(fac);
      par_lbls = lstPIs$piFEs$fe_pvs_lbls;
      names(lstRTMB_Params[[str]]) = par_lbls;
      names(lstRTMB_Map[[str]])    = par_lbls;
      rm(pars,npars,fac,par_lbls);
      dfrMM = lstRTMB_MMs[[str]] |> dplyr::distinct();
      rm(str);
    }#--if lstPIs$piFEs exists

    if (!is.null(lstPIs$piREs)){
      #--random effects: need to loop over each RE term (TODO!!)

    }#--lstPIs$piREs

    if (!is.null(lstPIs$piSMs)){
      #--smooths: need to loop over each smooth term (TODO!!)

    }#--lstPIs$piSMs

    lstMM[[par_id]] = dfrMM;
    rm(dfrMM,lstPIs,lstPVs);
  } #--par_ids loop

  #--expand model matrices in lstMM to largest coverage
  dfrMM = lstMM[[par_ids[1]]] |> dplyr::mutate(rownumber=dplyr::row_number());
  names(dfrMM)[names(dfrMM)=="rownumber"] = par_ids[1];
  if (length(lstMM)>1){
    for (i in 2:length(lstMM)){
      #--testing: i=2;
      cols0 = names(dfrMM);
      dfrMM1 = lstMM[[par_ids[i]]] |> dplyr::mutate(rownumber=dplyr::row_number());
      names(dfrMM1)[names(dfrMM1)=="rownumber"] = par_ids[i];
      cols1 = names(dfrMM1);
      dfrMM = dplyr::full_join(dfrMM,dfrMM1);
    }
  }


  #--extract parameter-level link functions
  vLinks = lst$dfrPEQs$link_fcn;      #--vector of link function names
  names(vLinks) = lst$dfrPEQs$par_id;

  lst = rtmbGMACS::namedList(
          lstRTMB_Params,   #--list with initial parameter values
          lstRTMB_TForms,   #--list with parameter value-level transform function names
          lstRTMB_LBs,      #--list with parameter value lower bounds
          lstRTMB_UBs,      #--list with parameter value upper bounds
          lstRTMB_Phases,   #--list with parameter value estimation phases
          lstRTMB_Mirrors,  #--list with parameter values mirroring
          lstRTMB_Xs,       #--list with (sparse) design matrices
          lstRTMB_MMs,      #--list with model matrices (as tibbles)
          lstRTMB_Map,      #--RTMB "map" list
          dfrMM,            #--fully-joined model_matrix with parameter indices (as tibble)
          vLinks            #--vector with parameter-level link function names

        );
  return(lst);
}
#--testing: lstRTMB = createRTMB_Params(lstCTL$lstAllometryInfo$lstParsInfo); #--input list has elements lstFcnInfo, lstParsInfo

