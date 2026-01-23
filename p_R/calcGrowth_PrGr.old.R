#'
#' @title Calculate growth
#' @description Function to calculate growth.
#' @param pGrA - mean post-molt size at zSclGrA
#' @param zGrA - pre-molt size at yielding pGrA as mean post-molt size
#' @param pGrB - mean post-molt size at zSclGrB
#' @param zGrB - pre-molt size at yielding pGrB as mean post-molt size
#' @param pGrBeta - gamma distribution scale parameter for post-molt variability
#' @param zBs - pre-molt sizes at which to calculate matrix
#' @return growth matrix
#'
#' @details The formula used is
#'
#' mnZs = pGrA\*exp(log(pGrB/pGrA)/log(zGrB/zGrA)\*log(zBs/zGrA));
#'
#' @examples
#' # example code
#' grM = grwPwrLaw1(33,25,150,125,1,seq(25,180,5));
#'
#' @md
#' @export
#'
grwPwrLaw1<-function(pGrA,zGrA,pGrB,zGrB,pGrBeta,zBs){
  mnZs = pGrA*exp(log(pGrB/pGrA)/log(zGrB/zGrA)*log(zBs/zGrA));
  return(grM);
}
#'
#' @title Calculate size-dependent growth
#' @description Function to calculate size-dependent growth.
#' @param pLnM - log-scale base mortality
#' @param pZ0 - reference size (fixed parameter)
#' @param zBs - pre-molt sizes at which to calculate growth
#' @return object with same dimensions as `z`.
#'
#' @details
#'
#' The formula for mean growth used is
#'
#' mnZs = exp(grA+grB*log(zBs));
#'
#' @examples
#' # example code
#' z = seq(25,100,5);
#' M = natMortZ(log(0.2),100,z);
#'
#' @md
#' @export
#'
grwPwrLaw2<-function(pGrA,zGrA,pGrB,zGrB,pGrBeta,zBs){
  mnZs = exp(grA+grB*log(zBs));
  return(grM);
}

#'
#' @title Calculate growth for all model categories across time
#' @description
#' Function to calculate growth for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Growth()])
#' @param params - RTMB parameters list with growth-specific elements
#' @param verbose - flag to print diagnostic info
#'
#' @return TODO: might want to return a list of a list of matrices
#'
#' @details Growth matrix for any given y_, s_ is technically an upper triangle
#' block-diagonal, with non-zero elements only for z_row <= z_column (i.e., post-molt size)
#' AND {r,x,m,p}_row == {r,x,m,p}_, followed by ??
#'
#' At start, {r,x,m,p,z} has probability of molting prM(r,x,m,p,z), which splits
#' n_{r,x,m,p,z} into molting (mn_{r,x,m,p,z}) and non-molting (nn_{r,x,m,p,z}) components.
#' If terminal molt depends on pre-molt size, it should be evaluated now on molting animals
#' (e.g. immature->mature for mn).
#'
#' The non-molting component should have p->max(p+1,p_max).
#'
#' The molting component undergoes growth as a block-diagonal with non-zero transitions possible only
#' for z_row <= z_column (i.e., post-molt size) AND {r,x,m,p}_row == {r,x,m,p}_column, followed
#' by p->post-molt age 0.
#'
#' If terminal molt depends on post-molt size, it would be evaluated now on molted crab
#' (e.g., immature-> mature).
#'
#' @export
#'
calGrowth<-function(dims,info,params,verbose=FALSE,loopIC_=TRUE){
  if (verbose) cat("Starting calcGrowth.\n")
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

    ###--calculate growth array----
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
    stop("unrecognized type option for growth:",info$option);
  }
  return(M);
}#--end of function
