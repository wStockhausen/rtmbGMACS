#'
#' @title Determine initial parameter values (revised)
#' @description Function (revised) to determine initial parameter values from the CTL list
#' and the model matrix list.
#' @param lstCTL - list created reading the CTL file
#' @param lstModMtx - list of model matrices created using [createParamsMapRev()]
#' @return list with elements `lst_pinfo`, `lst_priors`, `lst_map`, and `params`, which are in turn
#' named lists with information associated with each RTMB parameter.
#' @details
#' Additional details...
#'
#' @import dplyr
#' @export
#'
determineInitialValuesRev<-function(lstCTL,lstModMtx){
  #--determine RTMB parameters and initial values----
  dfrPEQs = lstCTL$dfrPEQs;
  dfrFEVs = lstCTL$dfrFEVs;
  dfrECVs = lstCTL$dfrECVs;
  dfrREVs = lstCTL$dfrREVs;
  #--create map for parameter phasing, mirroring
  ##--get priors info for model parameters
  ##--determine initial values for "actual" RTMB parameters
  ###--using svd inversion and link function to determine initial values
  ###--of "actual" RTMB parameters from input values for model parameters
  lst_pinfo  = list();
  lst_map    = list();
  lst_priors = list();
  params     = list();
  ##--process fixed effects info----
  if (!is.null(dfrFEVs)){
    for (par_id_ in unique(dfrFEVs$par_id)){
      ##--testing: par_id_="pLnA";
      dfrIVs   = dfrFEVs |> dplyr::filter(par_id==par_id_);               #--dataframe with model parameter info
      for (par_idx_ in unique(dfrIVs$par_idx)){
        #--testing: par_idx_ = unique(dfrIVs$par_idx)[1];
        rw  = dfrPEQs |> dplyr::filter(par_id==par_id_,par_idx==par_idx_);
        nm  = paste0("allom_",rw$par_key);#--RTMB parameter base name
        nmp = paste0(nm,"_FEs");
        dfrIVsp = dfrIVs |> dplyr::filter(par_idx_==par_idx) |> dplyr::mutate(par_nm=nmp);
        if ("tform" %in% names(dfrIVsp)){
          for (r in 1:nrow(dfrIVsp)){
            dfrIVsp$IV[r] = get(dfrIVsp$tform[r])(dfrIVsp$IV[r]);
          }
        }
        ##--create "actual" RTMB parameter vector corresponding to model parameter `par_id_`
        lst_pinfo[[nmp]]  = getParametersInitialValues(lst=lstModMtx[[nm]]$FEs,dfrIVsp,links=lstModMtx[[nm]]$links);
        params[[nmp]]     = lst_pinfo[[nmp]]$vP;
        lst_map[[nmp]]    = createParamsMapRev(dfrIVsp);       #--create RTMB "map" for model parameter vector
        lst_priors[[nmp]] = getParametersPriorInfo(dfrIVsp);#--get prior info for parameters
      } #--loop: par_idx_ in unique(dfrIVs$par_idx)
    } #--loop: par_id_ in names(dfrFEVs)
  } #--if: length(dfrFEVs)>0

  ##--process environmental covariates info----
  if (!is.null(dfrECVs)){
    for (par_id_ in unique(dfrECVs$par_id)){ #--loop over model parameter names
      ##--testing: par_id_ = unique(dfrECVs$par_id)[1];
      dfrIVs   = dfrECVs |> dplyr::filter(par_id==par_id_);  #--dataframe with model parameter info
      for (par_idx_ in unique(dfrIVs$par_idx)){ #--loop over "indices" associated with the model parameter
        #--testing: par_idx_ = unique(dfrIVs$par_idx)[1];
        rw  = dfrPEQs |> dplyr::filter(par_id==par_id_,par_idx==par_idx_);
        nm  = paste0("allom_",rw$par_key);#--RTMB parameter base name
        nmp = paste0(nm,"_ECs");
        dfrIVsp = dfrIVs |> dplyr::filter(par_idx_==par_idx) |> dplyr::mutate(par_nm=nmp);
        if ("tform" %in% names(dfrIVsp)){
          for (r in 1:nrow(dfrIVsp)){
            dfrIVsp$IV[r] = get(dfrIVsp$tform[r])(dfrIVsp$IV[r]);
          }
        }
        ##--create "actual" RTMB parameter vector corresponding to model parameter `par_id_`
        lst_pinfo[[nmp]]  = getParametersInitialValues(lst=lstModMtx[[nm]]$ECs,dfrIVsp=dfrIVsp,links=lstModMtx[[nm]]$links);
        params[[nmp]]     = (lst_pinfo[[nmp]])$vP;
        lst_map[[nmp]]    = createParamsMapRev(dfrIVsp);    #--create RTMB "map" for model parameter vector
        lst_priors[[nmp]] = getParametersPriorInfo(dfrIVsp);#--get prior info for parameters
      }#--par_idx_ in unique(dfrIVs$par_idx) loop
    }#--par_id_ in names(dfrECVs)
  }#--end processing dfrECVs

  ##--process random effects info  TODO!!----
  if (!is.null(dfrREVs)){
    ###-- TODO!!
  }

  return(list(lst_pinfo=lst_pinfo,lst_priors=lst_priors,lst_map=lst_map,params=params));
}

