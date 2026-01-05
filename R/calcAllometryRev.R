#'
#' @title Simple power law to calculate weight-at-size
#' @description Function to calculate weight-at-size using a simple power law.
#' @param z - size vector
#' @param pA - scaling parameter (scalar or same size vector as `z`)
#' @param pB - power parameter (scalar or same size vector as `z`)
#' @param pR - random effect value (scalar or same size vector as `z`)
#' @param link - link function for random effect ("none", "add" for additive effect, "mlt" for multiplicative effect)
#' @return vector with the predicted weights-at-size
#' @details Calculates weight `w` as
#' w = pA * z^pB
#'
#' @export
#'
pwrlaw1<-function(z,pA,pB,pR,link="none"){
  w = pA * (RTMB::AD(z)^pB);
  if (link!="none"){
    if (link=="add") {w = w + pR;} else
    if (link=="mlt") {w = w * pR;};
  }
  return(w);
}

#'
#' @title Alternative power law to calculate weight-at-size
#' @description Function to calculate weight-at-size using an alternative power law.
#' @param z - size vector
#' @param pLnA - log-scale scaling parameter (scalar or same size vector as `z`)
#' @param pB - power parameter (scalar or same size vector as `z`)
#' @param pZ0 -  reference size suc that w(pZ0) = exp(pLnA)
#' @param pR - random effect value (scalar or same size vector as `z`)
#' @param link - link function for random effect ("none", "add" for additive effect, "mlt" for multiplicative effect)
#' @return vector with the predicted weights-at-size
#' @details Calculates weight `w` as
#' w = exp(pLnA + pB*log(z/pZ0))
#'
#' @export
#'
pwrlaw2<-function(z,pLnA,pB,pZ0,pR,link="none"){
  w = exp(pLnA+pB*log(RTMB::AD(z)/pZ0));
  if (link!="none"){
    if (link=="add") {w = w + pR;} else
    if (link=="mlt") {w = w * pR;};
  }
  return(w);
}

#'
#' @title Predict allometric (weight-at-size) values given model-scale parameters (revised)
#' @description Function (revised) to predict allometric (weight-at-size) values given model-scale parameters.
#' @param dataDFR - dataframe with columns specifying allometric functions and allometric parameter indices
#' @param lst_modparams - list of model-scale parameters
#' @return dataDFR, with two columns added: "pars", a list column, and "prd", a numeric column. See Details.
#' @details The "pars" column of the returned dataframe is a list column with each row giving the values of the parameters
#' used to predict the weight for that row. The "ped" column is the predicted weight.
#'
#' @export
calcAllometryRev<-function(dataDFR,lst_modparams){
  dataDFR = dataDFR |> dplyr::mutate(pars=list(list()),
                                     prd=NA_real_);
  for (r in 1:nrow(dataDFR)){
    #testing: r = 2;
    rw = dataDFR[r,];
    if (rw$fcn_nm=="pwrLaw1"){
      #--evaluate pA
      lst_pA = rw$par_idxs[[1]]$pA;
      pA = RTMB::AD(0.0);
      nm = lst_pA$par_nm;
      if (!is.na(lst_pA$FEs)) pA = pA + (lst_modparams[[paste0(nm,"_FEs")]])[lst_pA$FEs];
      if (!is.na(lst_pA$ECs)) pA = pA + (lst_modparams[[paste0(nm,"_ECs")]])[lst_pA$ECs];
      if (!is.na(lst_pA$REs)) stop("not yet implemented!");
#      pA = link_inv(pA);#--TODO!!
      #--evaluate pB
      lst_pB = rw$par_idxs[[1]]$pB;
      pB = RTMB::AD(0.0);
      nm = lst_pB$par_nm;
      if (!is.na(lst_pB$FEs)) pB = pB + (lst_modparams[[paste0(nm,"_FEs")]])[lst_pB$FEs];
      if (!is.na(lst_pB$ECs)) pB = pB + (lst_modparams[[paste0(nm,"_ECs")]])[lst_pB$ECs];
      if (!is.na(lst_pB$REs)) stop("not yet implemented!");
#      pB1 = link_inv(pB1);#--TODO!!
      prd = pwrlaw1(RTMB::AD(as.numeric(rw$z)),pA=pA,pB=pB,pR=RTMB::AD(0.0));
      dataDFR$pars[r] = list(list(pA=pA,pB=pB,pR=0.0));#--no conversion necessary
    }  else if (rw$fcn_nm=="pwrLaw2"){
      #--evaluate pLnA
      lst_pLnA = rw$par_idxs[[1]]$pLnA;
      pLnA = RTMB::AD(0.0);
      nm = lst_pLnA$par_nm;
      if (!is.na(lst_pLnA$FEs)) pLnA = pLnA + (lst_modparams[[paste0(nm,"_FEs")]])[lst_pLnA$FEs];
      if (!is.na(lst_pLnA$ECs)) pLnA = pLnA + (lst_modparams[[paste0(nm,"_ECs")]])[lst_pLnA$ECs];
      if (!is.na(lst_pLnA$REs)) stop("not yet implemented!");
#      pLnA = link_inv(pLnA);#--TODO!!
      #--evaluate pB
      lst_pB = rw$par_idxs[[1]]$pB;
      pB = RTMB::AD(0.0);
      nm = lst_pB$par_nm;
      if (!is.na(lst_pB$FEs)) pB = pB + (lst_modparams[[paste0(nm,"_FEs")]])[lst_pB$FEs];
      if (!is.na(lst_pB$ECs)) pB = pB + (lst_modparams[[paste0(nm,"_ECs")]])[lst_pB$ECs];
      if (!is.na(lst_pB$REs)) stop("not yet implemented!");
#      pB1 = link_inv(pB1);#--TODO!!
      #--evaluate pZ0
      lst_pZ0 = rw$par_idxs[[1]]$pZ0;
      pZ0 = RTMB::AD(0.0);
      nm = lst_pZ0$par_nm;
      if (!is.na(lst_pZ0$FEs)) pZ0 = pZ0 + (lst_modparams[[paste0(nm,"_FEs")]])[lst_pZ0$FEs];
      if (!is.na(lst_pZ0$ECs)) pZ0 = pZ0 + (lst_modparams[[paste0(nm,"_ECs")]])[lst_pZ0$ECs];
      if (!is.na(lst_pZ0$REs)) stop("not yet implemented!");
#      pZ01 = link_inv(pZ01);#--TODO!!
      prd = pwrlaw2(RTMB::AD(as.numeric(rw$z)),pLnA=pLnA,pB=pB,pZ0=pZ0,pR=RTMB::AD(0.0));
      dataDFR$pars[r] = list(list(pLnA=pLnA,pB=pB,pR=0.0));#--no conversion necessary
    }  #--if: rw$fcn_id=="pwrLaw2"
    dataDFR$prd[r]  = prd;                    #--no conversion necessary
  }#--loop: r
  return(dataDFR);
}

