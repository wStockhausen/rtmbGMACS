#'
#' @title Calculate the size-independent probability of undergoing a molt
#' @description Function to size-independent calculate the probability of undergoing a molt.
#' @param prM - size-independent probability of undergoing a molt
#' @param zBs - pre-molt sizes at which to calculate vector
#' @return object same size as zBs
#'
#' @details The returned object is the same size/shape as `zBs`, with the `prM`
#' for all elements
#'
#' @examples
#' # example code
#' prM = prMolt_Constant(0.5,seq(25,100,5));
#'
#' @md
#' @export
#'
prMolt_Constant<-function(pM,zBs){
  prM = AD(array(pM,dim=length(zBs)));
  return(prM);
}
#'
#' @title Calculate a size-dependent, descending normal probability of undergoing a molt
#' @description Function to calculate a size-dependent growth, descending normal probability of undergoing a molt.
#' @param mdZ - size at which curve starts to descend from 1
#' @param wdZ - width of descent (standard deviation of normal curve)
#' @param zBs - pre-molt sizes at which to calculate vector
#' @param dZ - size bin width to use to set
#' @return object with same dimensions as `zBs`.
#'
#' @details
#'
#' The formula for mean growth used is
#'
#' $$prM(zBs < mdZ) = 1.0$$
#' $$prM(zBs \le mdZ) = exp(-0.5*((zBs-mdZ) \over wdZ)^2)$$
#'
#' @examples
#' # example code
#' z = seq(25,100,5);
#' prM = prMolt_DscNormal(55,30,z);
#'
#' @md
#' @export
#'
prMolt_DscNormal<-function(mdZ,wdZ,zBs){
  prM = AD(array(0,dim=length(zBs)));
  zBs = as.numeric(zBs);
  # cat("in prMolt_DscNormal\n")
  # cat("\t",zBs," \n");
  # cat("\t",exp(-0.5*((zBs-mdZ)/wdZ)^2)," \n");
  # cat("\t",squarewave_right(mdZ,zBs),"\n")
  prM = AD(exp(-0.5*((zBs-mdZ)/wdZ)^2)*squarewave_right(mdZ,zBs));
  # cat("\t",prM,"\n");
  return(prM);
}

#'
#' @title Calculate the probability of undergoing a molt for all model categories across time
#' @description
#' Function to calculate the probability of undergoing a molt for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Growth_PrMolt()])
#' @param params - RTMB parameters list with elements specific to the probability of undergoing a molt
#' @param verbose - flag to print diagnostic info
#'
#' @return TODO: might want to return a list of a list of matrices
#'
#' @details Application of the probability of molting vector, prM, to the
#' population vector during the growth season decomposes it into molting and
#' non-molting components.
#'
#' @import dplyr
#'
#' @md
#' @export
#'
calcGrowth_PrMolt<-function(dims,info,params,verbose=FALSE,loopIC_=TRUE){
  if (verbose) cat("Starting calcGrowth_PrMolt.\n")
  prM = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
  if (info$option=="data"){
    ##--"data" option----
    p = params$pPrMolt_FPs;#--vector of weights-at-size
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
        prM[iy_,is_,] = p[pidx];
      }#--is_ loop
    }#--iy_ loop
  } else if (tolower(info$option)=="function"){
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
        p = p + params$pPrMolt_MPs[dfrUCr$mpr_idx[1]];
      }
      if ((!is.null(dfrUCr$opr_idx))&&(!is.na(dfrUCr$opr_idx[1]))) {
        if (dfrUCr$op_type=="additive") {
          p = p + params$pPrMolt_OPs[dfrUCr$opr_idx[1]];
        } else {p = p * params$pPrMolt_OPs[dfrUCr$opr_idx[1]];}
      }
      if ((!is.null(dfrUCr$dpr_idx))&&(!is.na(dfrUCr$dpr_idx[1]))) {
        if (dfrUCr$dv_type=="additive") {
          p = p + params$pPrMolt_DPs[dfrUCr$dpr_idx[1]];
        } else {p = p * params$pPrMolt_DPs[dfrUCr$dpr_idx[1]];}
      }
      vals[rw] = p;
    }
    if (verbose){
      print(dfrUCs |> dplyr::mutate(vals_=vals));
    }

    ###--create index vector into `vals` using the index names from dfrUCs----
    idxVals = 1:length(vals);
    names(idxVals) = dfrUCs$idx;

    ###--calculate probability of molting vector----
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
            if (tolower(dfrIdxs$fcn)=="constant"){
              prM[iy_,is_,ic_] = prMolt_Constant(vals[idxVals[dfrIdxs$pPrM]],dfrIdxs$z);
            } else
            if (tolower(dfrIdxs$fcn)=="dscnormal"){
              prM[iy_,is_,ic_] = prMolt_DscNormal(vals[idxVals[dfrIdxs$pMdZ]],
                                                  vals[idxVals[dfrIdxs$pWdZ]],
                                                  dfrIdxs$z);
            }
          }#--ic_ loop
        }
        if (!loopIC_){
          dfrDims  = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |> dplyr::mutate(ic_=dplyr::row_number());
          dfrIdxsA = dfrDims |> dplyr::left_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z));
          dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)=="constant");
          if (nrow(dfrIdxs) > 0){
            ic_ = dfrIdxs$ic_;
            prM[iy_,is_,ic_] = prMolt_Constant(vals[idxVals[dfrIdxs$pPrM]],dfrIdxs$z);
          }
          dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)=="dscnormal");
          if (nrow(dfrIdxs) > 0){
            ic_ = dfrIdxs$ic_;
            prM[iy_,is_,ic_] = prMolt_DscNormal(vals[idxVals[dfrIdxs$pMdZ]],
                                                vals[idxVals[dfrIdxs$pWdZ]],
                                                dfrIdxs$z);
          }
        }
      }#--is_ loop
    }#--iy_ loop
  } else {
    stop("unrecognized type option for growth:",info$option);
  }
  return(prM);
}#--end of function
