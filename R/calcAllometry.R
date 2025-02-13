#'
#' @title Calculate allometry for all model categories across time
#' @description
#' Function to calculate allometryfor all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Allometry()])
#' @param params - RTMB parameters list with allometry-specific elements
#' @param verbose - flag to print diagnostic info
#'
#' @return list
#'
#' @export
#'
calcAllometry<-function(dims,info,params,verbose=FALSE){
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
                                              by = join_by(y, s, r, x, m, p, z));
        pidx = dfrIdxs$pidx;
        wAtZ[iy_,is_,] = p[pidx];
      }#--is_ loop
    }#--iy_ loop
  } else if (info$option=="function"){
    ##--"function" option----
    ###--calculate inputs to functions----
    ####--for each input parameter to a function, p = MP + OP + DP + ...
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
    dfrUCs$val = vals;#--just for testing (this is NOT AD-compatible)
    names(vals) = dfrUCs$idx;

    idxVals = 1:length(vals);
    names(idxVals) = dfrUCs$idx;

    pwrLaw1<-function(pA,pB,z){
      return(pA*(as.numeric(z)^pB));
    }
    pwrLaw2<-function(pLnA,pLnS,pB,pZ0,z){
      return(exp(pLnA-pLnS+pB*(as.numeric(z)-pZ0)));
    }

    for (iy_ in 1:dims$nYs){
      #--iy_ = 1;
      y_ = dims$y[iy_];
      for (is_ in 1:dims$nSs){
        #--is_= 1;
        s_ = dims$s[is_];
        for (ic_ in 1:dims$nCs){
          #--ic_ = 1;
          dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_))[ic_,];
          dfrIdxs = dfrDims |> dplyr::left_join(info$dfrHCs,by = join_by(y, s, r, x, m, p, z));
          if (dfrIdxs$fcn=="pwrLaw1"){
            wAtZ[iy_,is_,ic_] = pwrLaw1(vals[idxVals[dfrIdxs$pA]],vals[idxVals[dfrIdxs$pB]],dfrIdxs$z);
          } else
          if (dfrIdxs$fcn=="pwrLaw2"){
            wAtZ[iy_,is_,ic_] = pwrLaw2(vals[idxVals[dfrIdxs$pLnA]],vals[idxVals[dfrIdxs$pLnS]],
                                        vals[idxVals[dfrIdxs$pB]],vals[idxVals[dfrIdxs$pZ0]],dfrIdxs$z);
          }
        }#--ic_ loop
      }#--is_ loop
    }#--iy_ loop
  } else {
    stop("unrecognized type option for allometry:",info$option);
  }
  return(wAtZ);
}#--end of function
