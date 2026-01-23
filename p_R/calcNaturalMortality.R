#'
#' @title Calculate natural mortality
#' @description Function to calculate natural mortality.
#' @param pLnM - ln-scale natural mortality rate
#' @return mortality rate on arithmetic scale
#'
#' @details The returned object is the same dimensions and type as the
#' input `pLnM` object.The formula used is
#'
#' M = exp(pLnM)
#'
#' @examples
#' # example code
#' M = natMort(log(0.2));
#'
#' @md
#' @export
#'
natMort<-function(pLnM){
  return(exp(pLnM));
}
#'
#' @title Calculate size-dependent natural mortality
#' @description Function to calculate size-dependent natural mortality.
#' @param pLnM - log-scale base mortality
#' @param pZ0 - reference size (fixed parameter)
#' @param z - sizes at which to calculate mortality rates
#' @return object with same dimensions as `z`.
#'
#' @details The returned object is the same dimensions and type as the
#' input `z` object.If `pLnM` and `pZ0` are not scalars, they should be
#' conformable to the shape of `z`.
#'
#' The formula used is
#'
#' M(z) = exp{pLnM)*(pZ0/z)}
#'
#' @examples
#' # example code
#' z = seq(25,100,5);
#' M = natMortZ(log(0.2),100,z);
#'
#' @md
#' @export
#'
natMortZ<-function(pLnM,pZ0,z){
  return(exp(pLnM)*(pZ0/as.numeric(z)));
}

#'
#' @title Calculate natural mortality for all model categories across time
#' @description
#' Function to calculate natural mortality for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_NaturalMortality()])
#' @param params - RTMB parameters list with natural mortality-specific elements
#' @param verbose - flag to print diagnostic info
#'
#' @return 3-d array with dimensions `\[nYs, nSs,nCs\]`, where `nYs` is the number
#' of model years, `nSs` is the number of seasons/year, and `nCs` is the number of
#' population categories.
#'
#' @export
#'
calcNaturalMortality<-function(dims,info,params,verbose=FALSE,loopIC_=TRUE){
  if (verbose) cat("Starting calcNaturalMortality.\n")
  M = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
  if (info$option=="data"){
    ##--"data" option----
    p = params$pNM_FPs;#--vector of weights-at-size
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
        #   M[iy_,is_,ic_] = p[pidx];
        # }
        dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_));
        dfrIdxs = dfrDims |> dplyr::left_join(info$dfrDims2Pars,
                                              by = dplyr::join_by(y, s, r, x, m, p, z));
        pidx = dfrIdxs$pidx;
        M[iy_,is_,] = p[pidx];
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
        p = p + params$pNM_MPs[dfrUCr$mpr_idx[1]];
      }
      if (!is.na(dfrUCr$opr_idx[1])) {
        if (dfrUCr$op_type=="additive") {
          p = p + params$pNM_OPs[dfrUCr$opr_idx[1]];
        } else {p = p * params$pNM_OPs[dfrUCr$opr_idx[1]];}
      }
      if (!is.na(dfrUCr$dpr_idx[1])) {
        if (dfrUCr$dv_type=="additive") {
          p = p + params$pNM_DPs[dfrUCr$dpr_idx[1]];
        } else {p = p * params$pNM_DPs[dfrUCr$dpr_idx[1]];}
      }
      vals[rw] = p;
    }
    if (verbose){
      print(dfrUCs |> dplyr::mutate(vals_=vals));
    }

    ###--create index vector into `vals` using the index names from dfrUCs----
    idxVals = 1:length(vals);
    names(idxVals) = dfrUCs$idx;

    ###--calculate natural mortality array----
    ####--TODO: reorganize to speed up?? (see TMB email list discussions on assignment)
    ####--I think ic_ loop can be vectorized
    for (iy_ in 1:dims$nYs){
      #--iy_ = 1;
      y_ = dims$y[iy_];
      for (is_ in 1:dims$nSs){
        #--is_= 1;
        s_ = dims$s[is_];
        ####--two possible approaches to calculating M across pop categories
        ####--is one faster?
        if (loopIC_){
          for (ic_ in 1:dims$nCs){
            #--ic_ = 1;
            dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_))[ic_,];
            dfrIdxs = dfrDims |> dplyr::left_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z));
            if (tolower(dfrIdxs$fcn)=="natmort"){
              M[iy_,is_,ic_] = natMort(vals[idxVals[dfrIdxs$pLnM]]);
            } else
            if (tolower(dfrIdxs$fcn)=="natmortz"){
              M[iy_,is_,ic_] = natMortZ(vals[idxVals[dfrIdxs$pLnM]],
                                        vals[idxVals[dfrIdxs$pZ0]],
                                        dfrIdxs$z);
            }
          }#--ic_ loop
        }
        if (!loopIC_){
          dfrDims  = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |> dplyr::mutate(ic_=dplyr::row_number());
          dfrIdxsA = dfrDims |> dplyr::left_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z));
          dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)=="natmort");
          if (nrow(dfrIdxs) > 0){
            ic_ = dfrIdxs$ic_;
            M[iy_,is_,ic_] = natMort(vals[idxVals[dfrIdxs$pLnM]]);
          }
          dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)=="natmortz");
          if (nrow(dfrIdxs) > 0){
            ic_ = dfrIdxs$ic_;
            M[iy_,is_,ic_] = natMortZ(vals[idxVals[dfrIdxs$pLnM]],
                                      vals[idxVals[dfrIdxs$pZ0]],
                                      dfrIdxs$z);
          }
        }
      }#--is_ loop
    }#--iy_ loop
  } else {
    stop("unrecognized type option for natural mortality:",info$option);
  }
  return(M);
}#--end of function
