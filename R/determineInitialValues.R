determineInitialValues<-function(lstCTL,lstModMtx){
  #--determine RTMB parameters and initial values----
  dfrFcnParInfo = lstCTL$dfrFcnParInfo;
  lstFEVs = lstCTL$lstFEVs;
  lstECVs = lstCTL$lstECVs;
  lstREVs = lstCTL$lstREVs;
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
  if (length(lstFEVs)>0){
    for (par_id_ in names(lstFEVs)){
      ##--testing: par_id_ = names(lstFEVs)[1];
      dfrIVs   = lstFEVs[[par_id_]];               #--dataframe with model parameter info
      for (par_idx_ in unique(dfrIVs$par_idx)){
        #--testing: par_idx_ = unique(dfrIVs$par_idx)[1];
        rw = dfrFcnParInfo |> dplyr::filter(par_id==par_id_,par_idx==par_idx_);
        nm = paste0("allom_",rw$fcn_id,"_",rw$fcn_idx,"_",par_id_,"_",par_idx_);#--RTMB parameter base name
        dfrIVsp = dfrIVs |> dplyr::filter(par_idx_==par_idx) |> dplyr::mutate(par_nm=paste0(nm,"_FEs"),fcn_idx=rw$fcn_idx,fcn_id=rw$fcn_id,fcn_nm=rw$fcn_nm);
        ##--create "actual" RTMB parameter vector corresponding to model parameter `par_id_`
        nmp = paste0(nm,"_FEs");
        lst_pinfo[[nmp]]  = getParametersInitialValues(lst=lstModMtx[[nm]]$FEs,dfrIVsp,links=lstModMtx[[nm]]$links);
        params[[nmp]]     = lst_pinfo[[nmp]]$vP;
        lst_map[[nmp]]    = createParamsMap(dfrIVsp);       #--create RTMB "map" for model parameter vector
        lst_priors[[nmp]] = getParametersPriorInfo(dfrIVsp);#--get prior info for parameters
      } #--loop: par_idx_ in unique(dfrIVs$par_idx)
    } #--loop: par_id_ in names(lstFEVs)
  } #--if: length(lstFEVs)>0

  ##--process environmental covariates info----
  if (length(lstECVs)>0){
    for (par_id_ in names(lstECVs)){ #--loop over model parameter names
      ##--testing: par_id_ = names(lstECVs)[1];
      dfrIVs   = lstECVs[[par_id_]];  #--dataframe with model parameter info
      for (par_idx_ in unique(dfrIVs$par_idx)){ #--loop over "indices" associated with the model parameter
        #--testing: par_idx_ = unique(dfrIVs$par_idx)[1];
        rw = dfrFcnParInfo |> dplyr::filter(par_id==par_id_,par_idx==par_idx_);
        nm = paste0("allom_",rw$fcn_id,"_",rw$fcn_idx,"_",par_id_,"_",par_idx_);#--RTMB parameter base name
        dfrIVsp = dfrIVs |> dplyr::filter(par_idx_==par_idx) |> dplyr::mutate(par_nm=paste0(nm,"_ECs"),fcn_idx=rw$fcn_idx,fcn_id=rw$fcn_id,fcn_nm=rw$fcn_nm);
        ##--create "actual" RTMB parameter vector corresponding to model parameter `par_id_`
        nmp = paste0(nm,"_ECs");
        lst_pinfo[[nmp]]  = getParametersInitialValues(lst=lstModMtx[[nm]]$ECs,dfrIVsp=dfrIVsp,links=lstModMtx[[nm]]$links);
        params[[nmp]]     = (lst_pinfo[[nmp]])$vP;
        lst_map[[nmp]]    = createParamsMap(dfrIVsp);       #--create RTMB "map" for model parameter vector
        lst_priors[[nmp]] = getParametersPriorInfo(dfrIVsp);#--get prior info for parameters
      }#--par_idx_ in unique(dfrIVs$par_idx) loop
    }#--par_id_ in names(lstECVs)
  }#--end processing lstECVs

  ##--process random effects info  TODO!!----
  ###-- TODO!!

  return(list(lst_pinfo=lst_pinfo,lst_priors=lst_priors,lst_map=lst_map,params=params));
}

