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
  if (testing) {
    cat("starting objfun_EstimateAllometry.\n");
    cat("model parameters are:\n", paste0("\t",names(params),"\n"));
  }

  #--define allometric functions----
  pwrlaw1<-function(z,pA,pB,pR,link="none"){
    w = pA * (AD(z)^pB);
    if (link!="none"){
      if (link=="add") {w = w + pR;} else
      if (link=="mlt") {w = w * pR;};
    }
    return(w);
  }

  pwrlaw2<-function(z,pLnA,pB,pZ0,pR,link="none"){
    w = exp(pLnA+pB*ln(AD(z)/pZ0));
    if (link!="none"){
      if (link=="add") {w = w + pR;} else
      if (link=="mlt") {w = w * pR;};
    }
    return(w);
  }

  if (testing) {
    cat("RTMB params: \n")
    for (nm in names(params)){
      cat("\t",nm,": ",paste(params[[nm]],collapse=" "),"\n",sep="");
    }
  }

  #--predict weight-at-size
  ##--calculate current model parameter values
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
  RTMB::REPORT(lst_modparams);
  if (testing){
    cat("modparams = \n");
    for (nm in names(lst_modparams)){
      cat("\t",nm,":",lst_modparams[[nm]],"\n",sep=" ");
    }
  }

  ##--predict individual weight-at-size based on model parameter values
  nll = RTMB::AD(0.0)
  for (r in 1:nrow(inputs$lstAllom$dataDFR)){
    #testing: r = 2;
    rw = inputs$lstAllom$dataDFR[r,];
    if (rw$fcn_nm=="pwrLaw1"){
      #--evaluate pA
      lst_pA = rw$pA[[1]];
      pA = RTMB::AD(0.0);
      nm = paste0("allom_",rw$fcn_id,"_",rw$fcn_idx,"_pA_",lst_pA$par_idx);
      if (!is.na(lst_pA$FEs)) pA = pA + (lst_modparams[[paste0(nm,"_FEs")]])[lst_pA$FEs];
      if (!is.na(lst_pA$ECs)) pA = pA + (lst_modparams[[paste0(nm,"_ECs")]])[lst_pA$ECs];
      if (!is.na(lst_pA$REs)) stop("not yet implemented!");
#      pA = link_inv(pA);#--TODO!!
      #--evaluate pB1
      lst_pB1 = rw$pB1[[1]];
      pB1 = RTMB::AD(0.0);
      nm = paste0("allom_",rw$fcn_id,"_",rw$fcn_idx,"_pB1_",lst_pB1$par_idx);
      if (!is.na(lst_pB1$FEs)) pB1 = pB1 + (lst_modparams[[paste0(nm,"_FEs")]])[lst_pB1$FEs];
      if (!is.na(lst_pB1$ECs)) pB1 = pB1 + (lst_modparams[[paste0(nm,"_ECs")]])[lst_pB1$ECs];
      if (!is.na(lst_pB1$REs)) stop("not yet implemented!");
#      pB1 = link_inv(pB1);#--TODO!!
      prd = pwrlaw1(RTMB::AD(as.numeric(rw$z)),pA=pA,pB=pB1,pR=RTMB::AD(0.0));
      if (rw$obs>0) nll = nll - dlnorm(as.numeric(rw$obs),log(prd),1,log=TRUE);  #--should depend on "family"
      if (testing){
        #--keep as numeric values
        inputs$lstAllom$dataDFR$pars[r] = list(list(pA=pA,pB1=pB1));#--no conversion necessary
        inputs$lstAllom$dataDFR$prd[r]  = prd;                      #--no conversion necessary
      } else {
        #--convert to numeric values
        inputs$lstAllom$dataDFR$pars[r] = list(list(pA=RTMB:::getValues(pA),pB1=RTMB:::getValues(pB1)));
        inputs$lstAllom$dataDFR$prd[r]  = RTMB:::getValues(prd);
      }
    }#--if: rw$fcn_id=="pwrLaw1"
  }#--loop: r
  #--calculate numeric values of individual nll's
  obs = inputs$lstAllom$dataDFR$obs;#--
  prd = inputs$lstAllom$dataDFR$prd;
  idx = obs>0;
  inputs$lstAllom$dataDFR$nll[idx] = -stats::dlnorm(obs[idx],log(prd[idx]),1.0,log=TRUE);#--should depend on "family"

  REPORT(lst_modparams)
  dfrZWs = inputs$lstAllom$dataDFR |>
             dplyr::select(obs_id,y,x,m,p,z,obs,prd,nll) |>
             dplyr::mutate(z=as.numeric(z),obs=as.numeric(obs));
  REPORT(dfrZWs)

  if (testing) {
    cat("end objective function\n");
    return(list(nll=nll,dfrZWs=dfrZWs))
  }
  return(nll);
}


