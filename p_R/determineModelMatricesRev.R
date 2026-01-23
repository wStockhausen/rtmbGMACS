determineModelMatricesRev<-function(dfrPEQs){
  ##--determine model matrices for conversion from actual RTMB parameters to model parameters
  #dfrPEQs = lstCTL$dfrPEQs
  idxFEs = idxECs = idxREs = 0;
  lstModMtx = list();
  for (row in 1:nrow(dfrPEQs)){
    #--testing: row = 1;
    rw = dfrPEQs[row,];
    #--fixed effects
    lstModMtxFEs = calcModelMatrixFEs(txt=rw$feEQs,dfrMF=get(rw$par_frame),ctrs=get(rw$feContrasts));
    idxFEs = idxFEs + lstModMtxFEs$npars;
    lstModMtxFEs$idx_end = idxFEs;
    #--env. covariates
    lstModMtxECs = calcModelMatrixFEs(txt=rw$ecEQs,dfrMF=get(rw$par_frame),ctrs=NULL,verbose=TRUE);
    idxECs = idxECs + lstModMtxECs$npars;
    lstModMtxECs$idx_end = idxECs;
    #--random effects
    lstModMtxREs = calcModelMatrixREs(rw$reEQs,dfrMF=get(rw$par_frame),cov_type=rw$reCovType);
    idxREs = idxREs + lstModMtxREs$npars;
    lstModMtxREs$idx_end = idxREs;
    #--combine lists into one object
    nm=paste0("allom_",rw$par_key);
    lstModMtx[[nm]] = list(par_nm=nm,par_key=rw$par_key,par_id=rw$par_id,links=getLinkFcn(rw$link_fcn),
                           FEs=lstModMtxFEs,ECs=lstModMtxECs,REs=lstModMtxREs);
    rm(rw,lstModMtxFEs,lstModMtxECs,lstModMtxREs);
  }
  rm(row,nm);
  return(lstModMtx);
}
