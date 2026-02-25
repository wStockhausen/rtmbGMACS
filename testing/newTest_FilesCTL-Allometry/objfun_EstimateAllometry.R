#'
#' @title Calculate the objective function for a gmacs model
#' @description Function to calculate the objective function for a gmacs model.
#' @param params - named list of parameter values
#'
#' @return The objective function (negative log-likelihood) value given the data and parameters
#'
#' @details TODO!!
#' @examplesIf FALSE
#' # example code: TODO
#'
#'
#' @import RTMB
#'
#' @export
#' @md
#'
objfun_EstimateAllometry<-function(params){
  nll = RTMB::AD(0.0);
  #--testing: testing = TRUE;
  if (testing) {
    cat("starting objfun_EstimateAllometry.\n");
    cat("model parameters are:\n", paste0("\t",names(params),"\n"));
  }

  if (testing) {
    cat("RTMB params: \n")
    for (nm in names(params)){
      cat("\t",nm,": ",paste(params[[nm]],collapse=" "),"\n",sep="");
    }
    rm(nm);
  }

  #--predict weight-at-size
  ##--calculate current model parameter values from RTMB parameter values
  if (testing) cat("calculating current model parameter values.\n",sep="");
  lstAllomParamVals = list();
  nms_qpar = names(inputs$lstAllom$lstRTMB$lstRTMB_Params); #--"qualified" (by fe_, re_, or sm_) parameter names
  nms_upar = stringr::str_split_fixed(nms_qpar,"_",2)[,2];  #--"unqualified" parameter names
  if (testing) {
    cat("nms_upar =",nms_upar,"\n");
    cat("nms_qpar =",nms_qpar,"\n");
  }
  for (nm_upar in nms_upar){
    #--testing: nm_upar = nms_upar[1];
    if (testing) cat("\tcalculating values for ",nm_upar,".\n",sep="");
    ##--FEs:
    vs_fe = RTMB::AD(0.0);
    if (paste0("fe_",nm_upar) %in% nms_qpar){
      nm_qpar = paste0("fe_",nm_upar);
      if (testing) cat("\t\tcalculating FE values for ",nm_qpar,".\n",sep="");
      if (testing) cat("dim(X) =",dim(inputs$lstAllom$lstRTMB$lstRTMB_Xs[[nm_qpar]]),"length(p) =",length(params[[nm_qpar]]),"\n");
      vs_fe = as.vector(inputs$lstAllom$lstRTMB$lstRTMB_Xs[[nm_qpar]] %*% params[[nm_qpar]]);             #--FE value
      is_fe = as.vector(inputs$lstAllom$lstRTMB$lstRTMB_Xs[[nm_qpar]] %*% (1:length(params[[nm_qpar]]))); #--index of FE par level contributing
      if (testing) {
        cat("class(vs_fe) =",class(vs_fe),"\n");
        cat("length(vs_fe) =",length(vs_fe),"\n");
        cat(paste0("\tvalues for ",nm_qpar,":\n"),vs_fe,"\n",collapse=" ");
      }
    }
    ##--REs: TODO (do covariance and calc nll here, before any expansion of dimensions)
    vs_re = RTMB::AD(0.0);
    ##--SMs: TODO
    vs_sm = RTMB::AD(0.0);

    vs = vs_fe + vs_re + vs_sm;
    if (testing) cat(paste0("\tvalues for pre-link transform vs:\n"),vs,"\n",collapse=" ")
    ##--scale to function parameter value using (inverse?) link function?
    ###--`linkfcn` defined for function parameter, but `tform` defined for parameter levels
    ###--TODO: not convinced these should be consistent in all use cases (e.g., different units?), so leave as potentially different
    if (testing) cat("link function:",inputs$lstAllom$lstRTMB$vLinks[nm_upar],"\n");
    vs = tform_inv(inputs$lstAllom$lstRTMB$vLinks[nm_upar],vs);
    if (testing) cat(paste0("\tvalues for post-link transform vs:\n"),vs,"\n",collapse=" ")

    lstAllomParamVals[[nm_upar]] = vs;

  }#--loop over nm_upars

  if (testing){
    cat("lstAllomParamVals = \n");
    for (nm in names(lstAllomParamVals)){
      cat("\t",nm,":",lstAllomParamVals[[nm]],"\n",sep=" ");
    }
    rm(nm);
  }

  #--predict individual weight-at-size based on model parameter values
  dataDFR = calcAllometry(inputs$lstAllom$dataDFR,lst_modparams) |>
              dplyr::mutate(nll=NA_real_);
  ##--calculate nll for fit to data
  for (r in 1:nrow(dataDFR)){
    #testing: r = 2;
    rw = dataDFR[r,];
    if (rw$obs>0) {
      nllp = -dlnorm(as.numeric(rw$obs),log(rw$prd),1,log=TRUE);  #--should depend on "family"
      dataDFR$nll[r]  = nllp;                    #--no conversion necessary
      nll = nll + nllp;
    }
  }#--loop: r
  #
  # ##--output model predictions using ADREPORT to be able to obtain uncertainty intervals
  # for (nm in names(lst_modparams)){
  #   v = lst_modparams[[nm]];
  #   for (i in seq_along(v)) {
  #     eval(parse(text=paste0("modparam_",nm,"_",i,"=v[i];")));
  #     eval(parse(text=paste0("RTMB::ADREPORT(modparam_",nm,"_",i,")")));
  #   }
  # }
  #
  # ##--output fits to data using REPORT (no derivative information)
  # dfrZWs = dataDFR |>
  #            dplyr::select(obs_id,y,x,m,p,z,obs,prd,nll) |>
  #            dplyr::mutate(z=as.numeric(z),obs=as.numeric(obs));
  # RTMB::REPORT(dfrZWs)
  #
  # ##--calculate predicted values for factor/covariate combinations in inputs$lstAllom$dfrPrd
  # ###--output predictions using ADREPORT to be able to obtain uncertainty intervals
  # dfrPrd = calcAllometryRev(inputs$lstAllom$dfrPrd,lst_modparams);
  # prd_allom = dfrPrd$prd;
  # RTMB::ADREPORT(prd_allom);
  #
  # if (testing) {
  #   cat("end objective function\n");
  #   inputs$lstAllom$dataDFR = dataDFR;
  #   return(list(nll=nll,dfrZWs=dfrZWs,dfrPrd=dfrPrd))
  # }
  return(nll);
}


