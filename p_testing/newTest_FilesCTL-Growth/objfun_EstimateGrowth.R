#'
#' @title Calculate the objective function for a gmacs model
#' @description Function to calculate the objective function for a gmacs model.
#' @param params - named list of parameter values (p=1 for now)
#'
#' @return The objective function (negative log-likelihood) value given the data and parameters
#'
#' @details Input data and options should be organized in the calling environment as a list of lists
#' with the name "inputs". `inputs` should have the structure
#' \itemize{
#'   \item{dims - a list with model dimensions maps}
#'   \item{data - a named list with data-related entries}
#'   \item{priorsInfo - list with information on parameter-related priors}
#'   \item{testing - a logical to indicate the function is run in a testing environment}
#'   \item{testing - a logical to print intermediate results}
#' }
#'
#' Example
#' ```{r, eval=FALSE}
#' inputs = list(dims=NULL,data=NULL,priorInfo=NULL,testing=TRUE,testing=TRUE);
#' params = list(dummy=0.0);
#' obj = MakeADFun(objfun,params,silent=TRUE);
#' opt = nlminb(obj$par, obj$fn, obj$gr);
#' sdreport(obj);
#' ````
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

  pwrlaw2<-function(z,pLnA,pLnS,pB,pZ0,pR,link="none"){
    w = exp(pLnA+pLnS+pB*ln(AD(z)/pZ0));
    if (link!="none"){
      if (link=="add") {w = w + pR;} else
      if (link=="mlt") {w = w * pR;};
    }
    return(w);
  }

  if (testing) {
    cat("RTMB params: \n")
    for (nm in names(params)){
      cat("\t",nm,": ",params[[nm]],"\n",sep="");
    }
  }

  #--predict weight-at-size
  ##--calculate current model parameter values
  if (testing) cat("calculating current model parameter values.\n",sep="");
  lst_modparams = list();
  for (nm in names(inputs$lstAllom$lstModMtx)){
    #--testing: nm = names(inputs$lstAllom$lstModMtx)[1];
    if (testing) cat("\tmodel parameter values for ",nm,".\n",sep="");
    lstModMtx = inputs$lstAllom$lstModMtx[[nm]];
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
  for (r in 1:nrow(dfrZWppp)){
    #testing: r = 2;
    rw = dfrZWppp[r,];
    if (rw$`function`=="pwrLaw1"){
      #--evaluate pA
      lst_pA = rw$pA[[1]];
      pA = RTMB::AD(0.0);
      nm = paste0("allom_pA-",lst_pA$par_idx);
#      link_inv = inputs$lstAllom$lstModMtx[[nm]]$links$link_inv;#TODO: is this an AD function at this point?
      # if (!is.na(lst_pA$FEs)) pA = pA + (inputs$lstAllom$lstModMtx[[nm]]$FEs$mtxRedMM %*% params[[paste0(nm,"_FEs")]])[lst_pA$FEs];
      # if (!is.na(lst_pA$ECs)) pA = pA + (inputs$lstAllom$lstModMtx[[nm]]$ECs$mtxRedMM %*% params[[paste0(nm,"_ECs")]])[lst_pA$ECs];
      if (!is.na(lst_pA$FEs)) pA = pA + (lst_modparams[[paste0(nm,"_FEs")]])[lst_pA$FEs];
      if (!is.na(lst_pA$ECs)) pA = pA + (lst_modparams[[paste0(nm,"_ECs")]])[lst_pA$ECs];
      if (!is.na(lst_pA$REs)) stop("not yet implemented!");
#      pA = link_inv(pA);
      #--evaluate pB1
      lst_pB1 = rw$pB1[[1]];
      pB1 = RTMB::AD(0.0);
      nm = paste0("allom_pB1-",lst_pB1$par_idx);
      # if (!is.na(lst_pB1$FEs)) pB1 = pB1 + (inputs$lstAllom$lstModMtx[[nm]]$FEs$mtxRedMM %*% params[[paste0(nm,"_FEs")]])[lst_pB1$FEs];
      # if (!is.na(lst_pB1$ECs)) pB1 = pB1 + (inputs$lstAllom$lstModMtx[[nm]]$ECs$mtxRedMM %*% params[[paste0(nm,"_ECs")]])[lst_pB1$ECs];
      if (!is.na(lst_pB1$FEs)) pB1 = pB1 + (lst_modparams[[paste0(nm,"_FEs")]])[lst_pB1$FEs];
      if (!is.na(lst_pB1$ECs)) pB1 = pB1 + (lst_modparams[[paste0(nm,"_ECs")]])[lst_pB1$ECs];
      if (!is.na(lst_pB1$REs)) stop("not yet implemented!");
#      pB1 = link_inv(pB1);
      prd = pwrlaw1(RTMB::AD(as.numeric(rw$z)),pA=pA,pB=pB1,pR=RTMB::AD(0.0));
      if (rw$obs>0) nll = nll - dlnorm(as.numeric(rw$obs),log(prd),1,log=TRUE);
      dfrZWppp$pars[r] = list(list(pA=as.numeric(pA),pB1=as.numeric(pB1)));
      dfrZWppp$prd[r]  = as.numeric(prd);
    }#--if: rw$`function`=="pwrLaw1"
  }#--loop: r
  obs = as.numeric(dfrZWppp$obs);
  idx = obs>0;
  dfrZWppp$nll[idx] = -dlnorm(obs[idx],log(dfrZWppp$prd[idx]),1,log=TRUE);

  REPORT(lst_modparams)
  dfrZWs = dfrZWppp |> dplyr::select(obs_id,y,x,m,p,z,obs,prd,nll) |>
             dplyr::mutate(z=as.numeric(z),obs=as.numeric(obs));
  REPORT(dfrZWs)

  if (testing) cat("end objective function\n");
  return(nll);
}


