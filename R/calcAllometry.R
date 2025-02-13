#'
#' @title Calculate weight-at-size
#' @description Function to calculate weight-at-size.
#' @param pA - weight (kg) when z = 1 mm
#' @param pB - exponent on size
#' @param z - sizes at which to calculate weights (mm-scale)
#' @return weight in kg
#'
#' @details The returned object is the same dimensions and type as the
#' input `z` object.The formula used is
#'
#' W(z) = pA * z^{pB}
#'
#' such that W(z=1mm) = pA in kg;
#'
#' @examples
#' # example code
#' z = seq(25,100,5);
#' wAtZ = allom_pwrLaw1(0.00027/1000,3.022134,z)
#'
#' @md
#' @export
#'
allom_pwrLaw1<-function(pA,pB,z){
  return(pA*(as.numeric(z)^pB));
}
#'
#' @title Calculate weight-at-size
#' @description Function to calculate weight-at-size.
#' @param pLnA - log-scale weight when z = pZ0
#' @param pLnS - log-scale units conversion to kg (fixed parameter)
#' @param pB - exponent on relative size
#' @param pZ0 - reference size (fixed parameter)
#' @param z - sizes at which to calculate weights
#' @return weight in kg
#'
#' @details The returned object is the same dimensions and type as the
#' input `z` object.The formula used is
#'
#' W(z) = exp{pLnA-pLnS + pB*log(z/pZ0)}
#'
#' such that W(z=pZ0) = exp(pLnA)/exp(pLnS);
#'
#' @examples
#' # example code
#' z = seq(25,100,5);
#' wAtZ = allom_pwrLaw2(20,1000,2.0,50.0,z)
#'
#' @md
#' @export
#'
allom_pwrLaw2<-function(pLnA,pLnS,pB,pZ0,z){
  return(exp(pLnA-pLnS+pB*log(as.numeric(z)/pZ0)));
}

#'
#' @title Calculate allometry for all model categories across time
#' @description
#' Function to calculate allometryfor all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Allometry()])
#' @param params - RTMB parameters list with allometry-specific elements
#' @param verbose - flag to print diagnostic info
#'
#' @return 3-d array with dimensions `\[nYs, nSs,nCs\]`, where `nYs` is the number
#' of model years, `nSs` is the number of seasons/year, and `nCs` is the number of
#' population categories.
#'
#' @export
#'
calcAllometry<-function(dims,info,params,verbose=FALSE,loopIC_=TRUE){
  if (verbose) cat("Starting calcAllometry.\n")
  wAtZ = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
  if (info$option=="data"){
    ##--"data" option----
    p = params$pAllom_FPs;#--vector of weights-at-size
    #--need to expand to p to all years, seasons, and population categories
    for (iy_ in 1:dims$nYs){
      #--iy_ = 1;
      y_ = dims$y[iy_];
      for (is_ in 1:dims$nSs){
        #--is_= 1;
        s_ = dims$s[is_];
        # for (ic_ in 1:dims$nCs){
        #   #--ic_ = 1;
        #   dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_))[ic_,];
        #   dfrIdxs = dfrDims |> dplyr::left_join(info$dfrDims2Pars);
        #   pidx = dfrIdxs$pidx[1];
        #   wAtZ[iy_,is_,ic_] = p[pidx];
        # }
        dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_));
        dfrIdxs = dfrDims |> dplyr::left_join(info$dfrDims2Pars,
                                              by = dplyr::join_by(y, s, r, x, m, p, z));
        pidx = dfrIdxs$pidx;
        wAtZ[iy_,is_,] = p[pidx];
      }#--is_ loop
    }#--iy_ loop
  } else if (info$option=="function"){
    ##--"function" option----
    ###--calculate inputs to functions----
    ####--for each input parameter to a function, p = MP + OP + DP + ...
    ####--calculate values only for unique combinations of MP, OP, DP, etc.
    dfrUCs = info$dfrUniqCmbs;
    nRWs = nrow(dfrUCs);
    vals  = AD(array(0,nRWs)); #
    for (rw in 1:nrow(dfrUCs)){
      #--testing: rw = 1;
      dfrUCr = dfrUCs[rw,];
      p = AD(0);
      if (!is.na(dfrUCr$mpr_idx[1])) {
        p = p + params$pAllom_MPs[dfrUCr$mpr_idx[1]];
      }
      if (!is.na(dfrUCr$opr_idx[1])) {
        if (dfrUCr$op_type=="additive") {
          p = p + params$pAllom_OPs[dfrUCr$opr_idx[1]];
        } else {p = p * params$pAllom_OPs[dfrUCr$opr_idx[1]];}
      }
      if (!is.na(dfrUCr$dpr_idx[1])) {
        if (dfrUCr$dv_type=="additive") {
          p = p + params$pAllom_DPs[dfrUCr$dpr_idx[1]];
        } else {p = p * params$pAllom_DPs[dfrUCr$dpr_idx[1]];}
      }
      vals[rw] = p;
    }

    ###--create index vector into `vals` using the index names from dfrUCs----
    idxVals = 1:length(vals);
    names(idxVals) = dfrUCs$idx;

    ###--calculate weight-at-size array----
    ####--TODO: reorganize to speed up?? (see TMB email list discussions on assignment)
    ####--I think ic_ loop can be vectorized
    for (iy_ in 1:dims$nYs){
      #--iy_ = 1;
      y_ = dims$y[iy_];
      for (is_ in 1:dims$nSs){
        #--is_= 1;
        s_ = dims$s[is_];
        ####--two possible approaches to calculating wAtZ across pop categories
        ####--is one faster?
        if (loopIC_){
          for (ic_ in 1:dims$nCs){
            #--ic_ = 1;
            dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_))[ic_,];
            dfrIdxs = dfrDims |> dplyr::left_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z));
            if (dfrIdxs$fcn=="pwrLaw1"){
              wAtZ[iy_,is_,ic_] = allom_pwrLaw1(vals[idxVals[dfrIdxs$pA]],vals[idxVals[dfrIdxs$pB]],dfrIdxs$z);
            } else
            if (dfrIdxs$fcn=="pwrLaw2"){
              wAtZ[iy_,is_,ic_] = allom_pwrLaw2(vals[idxVals[dfrIdxs$pLnA]],vals[idxVals[dfrIdxs$pLnS]],
                                          vals[idxVals[dfrIdxs$pB]],vals[idxVals[dfrIdxs$pZ0]],dfrIdxs$z);
            }
          }#--ic_ loop
        }
        if (!loopIC_){
          dfrDims  = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |> dplyr::mutate(ic_=dplyr::row_number());
          dfrIdxsA = dfrDims |> dplyr::left_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z));
          dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)=="pwrlaw1");
          if (nrow(dfrIdxs) > 0){
            ic_ = dfrIdxs$ic_;
            wAtZ[iy_,is_,ic_] = allom_pwrLaw1(vals[idxVals[dfrIdxs$pA]],
                                              vals[idxVals[dfrIdxs$pB]],
                                              dfrIdxs$z);
          }
          dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)=="pwrlaw2");
          if (nrow(dfrIdxs) > 0){
            ic_ = dfrIdxs$ic_;
            wAtZ[iy_,is_,ic_] = allom_pwrLaw2(vals[idxVals[dfrIdxs$pLnA]],
                                             vals[idxVals[dfrIdxs$pLnS]],
                                             vals[idxVals[dfrIdxs$pB]],
                                             vals[idxVals[dfrIdxs$pZ0]],
                                             dfrIdxs$z);
          }
        }
      }#--is_ loop
    }#--iy_ loop
  } else {
    stop("unrecognized type option for allometry:",info$option);
  }
  return(wAtZ);
}#--end of function
