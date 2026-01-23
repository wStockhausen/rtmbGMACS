determineModelMatrices<-function(dfrFcnParInfo){
  ##--determine model matrices for conversion from actual RTMB parameters to model parameters
  #dfrFcnParInfo = dfrFcns |> dplyr::full_join(dfrPEQs,by="fcn_idx");
  idxFEs = idxECs = idxREs = 0;
  lstModMtx = list();
  for (par_idx_ in dfrFcnParInfo$par_idx){
    #--testing: par_idx_ = 1;
    rw = dfrFcnParInfo |> dplyr::filter(par_idx==par_idx_);
    #--fixed effects
    lstModMtxFEs = calcModelMatrixFEs(txt=rw$feEQs,dfrMF=get(rw$frame),ctrs=get(rw$feContrasts));
    idxFEs = idxFEs + lstModMtxFEs$npars;
    lstModMtxFEs$idx_end = idxFEs;
    #--env. covariates
    lstModMtxECs = calcModelMatrixFEs(txt=rw$ecEQs,dfrMF=get(rw$frame),ctrs=NULL,verbose=TRUE);
    idxECs = idxECs + lstModMtxECs$npars;
    lstModMtxECs$idx_end = idxECs;
    #--random effects
    lstModMtxREs = calcModelMatrixREs(rw$reEQs,dfrMF=get(rw$frame),cov_type=rw$reCovType);
    idxREs = idxREs + lstModMtxREs$npars;
    lstModMtxREs$idx_end = idxREs;
    #--combine lists into one object
    nm=paste0("allom_",rw$fcn_id,"_",rw$fcn_idx,"_",rw$par_id,"_",rw$par_idx);
    lstModMtx[[nm]] = list(par_nm=nm,fcn_idx=rw$fcn_idx,fcn_id=rw$fcn_id,fcn_nm=rw$fcn_nm,
                           par_idx=rw$par_idx,par_id=rw$par_id,links=getLinkFcn(rw$link_fcn),
                           FEs=lstModMtxFEs,ECs=lstModMtxECs,REs=lstModMtxREs);
    rm(rw,lstModMtxFEs,lstModMtxECs,lstModMtxREs);
  }
  rm(par_idx_,nm);
  return(lstModMtx);
}
