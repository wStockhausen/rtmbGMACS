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
objfun_EstimateAllometryRev<-function(params){
  #--testing: testing = TRUE;
  if (testing) {
    cat("starting objfun_EstimateAllometryRev.\n");
    cat("model parameters are:\n", paste0("\t",names(params),"\n"));
  }

  if (testing) {
    cat("RTMB params: \n")
    for (nm in names(params)){
      cat("\t",nm,": ",paste(params[[nm]],collapse=" "),"\n",sep="");
    }
  }

  #--predict weight-at-size
  ##--calculate current model parameter values from RTMB parameter values
  if (testing) cat("calculating current model parameter values.\n",sep="");
  lst_modparams = list();
  for (nm in names(inputs$lstModMtx)){
    #--testing: nm = names(inputs$stModMtx)[1];
    if (testing) cat("\tmodel parameter values for ",nm,".\n",sep="");
    lstModMtx = inputs$lstModMtx[[nm]];
    ##--FEs
    if (lstModMtx$FEs$npars>0){
      rtmbPars = params[[paste0(nm,"_FEs")]];
      T = as.matrix(lstModMtx$FEs$mtxRedMM0);
      modParFEs = as.vector(T %*% rtmbPars);
    } else {modParFEs=RTMB::AD(0.0);}
    lst_modparams[[paste0(nm,"_FEs")]] = modParFEs;
    if (testing) cat("\t\tmodParFEs =",modParFEs,"\n",sep=" ");
    ##--ECs
    if (lstModMtx$ECs$npars>0){
      rtmbPars = params[[paste0(nm,"_ECs")]];
      T = as.matrix(lstModMtx$ECs$mtxRedMM0);
      modParECs = as.vector(T %*% rtmbPars);
    } else {modParECs=RTMB::AD(0.0);}
    lst_modparams[[paste0(nm,"_ECs")]] = modParECs;
    if (testing) cat("\t\tmodParECs =",modParECs,"\n",sep=" ");
    ##--REs
    if (lstModMtx$REs$npars>0){
      rtmbPars = params[[paste0(nm,"_REs")]];
      T = as.matrix(lstModMtx$REs$mtxRedMM0);
      modParREs = as.vector(T %*% rtmbPars);
    } else {modParREs=RTMB::AD(0.0);}
    lst_modparams[[paste0(nm,"_REs")]] = modParREs;
    if (testing) cat("\t\tmodParREs =",modParREs,"\n",sep=" ");
    if (testing) cat(length(lst_modparams),object.size(lst_modparams),"\n\n")
#    rm(lstModMtx,modParFEs,modParECs,modParREs);
  }

  if (testing){
    cat("modparams = \n");
    for (nm in names(lst_modparams)){
      cat("\t",nm,":",lst_modparams[[nm]],"\n",sep=" ");
    }
  }

  ##--predict individual weight-at-size based on model parameter values
  nll = RTMB::AD(0.0);
  dataDFR = calcAllometryRev(inputs$lstAllom$dataDFR,lst_modparams) |>
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

  ##--output model predictions using ADREPORT to be able to obtain uncertainty intervals
  for (nm in names(lst_modparams)){
    v = lst_modparams[[nm]];
    for (i in seq_along(v)) {
      eval(parse(text=paste0("modparam_",nm,"_",i,"=v[i];")));
      eval(parse(text=paste0("RTMB::ADREPORT(modparam_",nm,"_",i,")")));
    }
  }

  ##--output fits to data using REPORT (no derivative information)
  dfrZWs = dataDFR |>
             dplyr::select(obs_id,y,x,m,p,z,obs,prd,nll) |>
             dplyr::mutate(z=as.numeric(z),obs=as.numeric(obs));
  RTMB::REPORT(dfrZWs)

  ##--calculate predicted values for factor/covariate combinations in inputs$lstAllom$dfrPrd
  ###--output predictions using ADREPORT to be able to obtain uncertainty intervals
  dfrPrd = calcAllometryRev(inputs$lstAllom$dfrPrd,lst_modparams);
  prd_allom = dfrPrd$prd;
  RTMB::ADREPORT(prd_allom);

  if (testing) {
    cat("end objective function\n");
    inputs$lstAllom$dataDFR = dataDFR;
    return(list(nll=nll,dfrZWs=dfrZWs,dfrPrd=dfrPrd))
  }
  return(nll);
}


