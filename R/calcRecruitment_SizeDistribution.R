#' #'
#' #' @title Calculate a gamma-distributed size distribution at recruitment
#' #' @description Function to calculate a gamma-distributed size distribution at recruitment.
#' #' @param mnZ - mean size at recruitment
#' #' @param wdZ - width (std. dev of gamma distribution)
#' #' @param zMn - minimum size with non-zero probability of recruitment
#' #' @param zMx - maximum size with non-zero probability of recruitment
#' #' @param zBs - sizes at which to calculate vector
#' #' @param dZ - size bin width to use to set
#' #' @return object with same dimensions as `zBs`.
#' #'
#' #' @details
#' #'
#' #' The formula used is
#' #'
#' #' $$prZ(zBs < zMn or zBs > zMx)   = 0.0$$
#' #' otherwise, for mnZ <=zBs<=zMx,
#' #' $$prZmn = pgamma(mnZ-dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ)$$
#' #' $$prZmx = pgamma(zMx+dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ)$$
#' #' $$prZ(zBs) = pgamma(zBs+dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ) -
#' #'              pgamma(zBs-dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ)$$
#' #' $$prZ(zBs) = prZ(zBs \le zMx)/(prZMx-prZmn)$$
#' #'
#' #' @examples
#' #' # example code
#' #' zBs = seq(25,100,5);
#' #' dZ  = z[2]-z[1];
#' #' prZ = prRecZ1(35,20,55,25,80,zBs,dZ);
#' #'
#' #' @md
#' #' @export
#' #'
#' prRecZ1<-function(mnZ,wdZ,zMn,zMx,zBs,dZ){
#'   prZ = AD(array(0,dim=length(zBs)));
#'   zBs = as.numeric(zBs);
#'   # cat("dZ:",dZ,"\n");
#'   # cat("ZBs\n");
#'   # print(zBs);
#'   # cat("mnZ\n");
#'   # print(mnZ);
#'   # cat("wdZ\n");
#'   # print(wdZ);
#'   # cat("zMn\n");
#'   # print(zMn);
#'   # cat("zMx\n");
#'   # print(zMx);
#'   #--calculate squarewave to "box" size range
#'   sqw = squarewave(zMn,zMx,zBs,dZ)
#'   # cat("sqw\n");
#'   # print(sqw);
#'   # cat("1-sqw\n");
#'   # print(1-sqw);
#'   shp = (mnZ^2)/(wdZ^2);
#'   # cat("shp\n");
#'   # print(shp);
#'   scl = (wdZ^2)/mnZ;
#'   # cat("scl\n");
#'   # print(scl);
#'   prZmn = pgamma(zMn-dZ/2,shape=shp,scale=scl);
#'   prZmx = pgamma(zMx+dZ/2,shape=shp,scale=scl);
#'   prZ   = pgamma(zBs+dZ/2,shape=shp,scale=scl) -
#'             pgamma(zBs-dZ/2,shape=shp,scale=scl);
#'   # cat("prZ before scaling\n")
#'   # print(prZ);
#'   # print(prZmx- prZmn);
#'   #--normalize and apply squarewave window so sum  across zBs is 1 for each population category
#'   prZ = sqw*prZ/(prZmx- prZmn);
#'   # cat("prZ after scaling\n")
#'   # print(prZ);
#'   # cat("sum(prZ):",sum(prZ),"\n");
#'   return(prZ);
#' }

#'
#' @title Calculate a gamma-distributed size distribution at recruitment
#' @description Function to calculate a gamma-distributed size distribution at recruitment.
#' @param mnZ - mean size at recruitment
#' @param wdZ - width (std. dev of gamma distribution)
#' @param zMn - minimum size with non-zero probability of recruitment (see Details)
#' @param zMx - maximum size with non-zero probability of recruitment (see Details)
#' @param zBs - sizes at which to calculate vector
#' @param dZ - size bin width to use to set
#' @return object with same dimensions as `zBs`.
#'
#' @details In an (RTMB) AD context, zMn and zMx must be fixed parameters identified in
#' the `map` list when ADMakeFun'ing an objective function that uses this function.
#'
#' The formula used is
#'
#' $$prZ(zBs < zMn or zBs > zMx)   = 0.0$$
#' otherwise, for mnZ <=zBs<=zMx,
#' $$prZmn = pgamma(mnZ-dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ)$$
#' $$prZmx = pgamma(zMx+dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ)$$
#' $$prZ(zBs) = pgamma(zBs+dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ) -
#'              pgamma(zBs-dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ)$$
#' $$prZ(zBs) = prZ(zBs)/(prZMx-prZmn)$$
#'
#' @examples
#' # example R code
#' zBs = seq(25,100,5);
#' dZ  = z[2]-z[1];
#' prZ = prRecZ1(35,20,55,25,80,zBs,dZ);
#'
#' @md
#' @export
#'
prRecZ1<-function(mnZ,wdZ,zMn,zMx,zBs,dZ){
  prZ = AD(array(0,dim=length(zBs)));
  if (RTMB:::ad_context()){
    nzMn = RTMB:::getValues(zMn);#--zMn is fixed, so this should not be a problem
    nzMx = RTMB:::getValues(zMx);#--zMx is fixed, so this should not be a problem
  } else {
    nzMn = zMn;
    nzMx = zMx;
  }
  idx = which((nzMn<=zBs)&(zBs<=nzMx));
  shp = (mnZ^2)/(wdZ^2);
  scl = (wdZ^2)/mnZ;
  prZmn = pgamma(zMn-dZ/2,shape=shp,scale=scl);
  prZmx = pgamma(zMx+dZ/2,shape=shp,scale=scl);
  prZ[idx] = pgamma(zBs[idx]+dZ/2,shape=shp,scale=scl) -
                 pgamma(zBs[idx]-dZ/2,shape=shp,scale=scl);
  # cat("prZ before scaling\n")
  # print(prZ);
  # print(prZmx- prZmn);
  #--normalize so sum  across zBs is 1 for each population category
  prZ[idx] = prZ[idx]/(prZmx- prZmn);
  # cat("prZ after scaling\n")
  # print(prZ);
  # cat("sum(prZ):",sum(prZ),"\n");
  return(prZ);
}

#'
#' @title Calculate the size distribution at recruitment for all model categories across time
#' @description
#' Function to calculate the size distribution at recruitment for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Recruitment_SizeDistribution()])
#' @param params - RTMB parameters list with elements specific to the probability of undergoing a molt
#' @param verbose - flag to print diagnostic info
#'
#' @return TODO: might want to return a list of a list of matrices
#'
#' @details Combination of the size distribution at recruitment vector, prZ, with other aspects of
#' recruitment yields recruitment by year, season, population category (excluding size) and size.
#'
#' @import dplyr
#'
#' @md
#' @export
#'
calcRecruitment_SizeDistribution<-function(dims,info,params,verbose=FALSE,loopIC_=FALSE){
  if (verbose) cat("Starting calcRecruitment_SizeDistribution.\n")
  prZ = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
  if (info$option=="data"){
    ##--"data" option----
    p = params$pRecZ_FPs;#--vector of recruitment size distribution values
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
        prZ[iy_,is_,] = p[pidx];
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
        p = p + params$pRecZ_MPs[dfrUCr$mpr_idx[1]];
      }
      if (any(names(dfrUCr)=="opr_idx"))
        if(!is.na(dfrUCr$opr_idx[1])) {
          if (dfrUCr$op_type=="additive") {
            p = p + params$pRecZ_OPs[dfrUCr$opr_idx[1]];
          } else {p = p * params$pRecZ_OPs[dfrUCr$opr_idx[1]];}
        }
      if (any(names(dfrUCr)=="dpr_idx"))
        if (!is.na(dfrUCr$dpr_idx[1])) {
          if (dfrUCr$dv_type=="additive") {
            p = p + params$pRecZ_DPs[dfrUCr$dpr_idx[1]];
          } else {p = p * params$pRecZ_DPs[dfrUCr$dpr_idx[1]];}
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
    dZ = unname(dims$zb[2]-dims$zb[1]);#--bin size
    for (iy_ in 1:dims$nYs){
      #--iy_ = 1;
      y_ = dims$y[iy_];
      for (is_ in 1:dims$nSs){
        #--is_= 1;
        s_ = dims$s[is_];
        ####--two possible approaches to calculating M across pop categories
        ####--is one faster?
        if (loopIC_){
          dmsC = dims$dmsYSC |> dplyr::filter(y==y_,s==s_);
          for (ic_ in 1:dims$nCs){
            #--ic_ = 1;
            dfrDims = dmsC[ic_,];
            dfrIdxs = dfrDims |> dplyr::inner_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z));
            if (nrow(dfrIdxs)>0){
              if (tolower(dfrIdxs$fcn)==tolower("prRecZ1")){
                prZ[iy_,is_,ic_] = prRecZ1(vals[idxVals[dfrIdxs$pMnZ]],
                                           vals[idxVals[dfrIdxs$pWdZ]],
                                           vals[idxVals[dfrIdxs$pZmn]],
                                           vals[idxVals[dfrIdxs$pZmx]],
                                           dfrIdxs$z,
                                           dZ);
              }
            }
          }#--ic_ loop
        }
        if (!loopIC_){
          dfrDims  = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |> dplyr::mutate(ic_=dplyr::row_number());
          dfrIdxsA = dfrDims |> dplyr::inner_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z));
          dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)==tolower("prRecZ1"));
          if (nrow(dfrIdxs) > 0){
            ic_ = dfrIdxs$ic_;
            prZ[iy_,is_,ic_] = prRecZ1(vals[idxVals[dfrIdxs$pMnZ]],
                                       vals[idxVals[dfrIdxs$pWdZ]],
                                       vals[idxVals[dfrIdxs$pZmn]],
                                       vals[idxVals[dfrIdxs$pZmx]],
                                       dfrIdxs$z,
                                       dZ);
          }
          #--other functions?
        }
      }#--is_ loop
    }#--iy_ loop
  } else {
    stop("unrecognized type option for calcRecruitment_SizeDistribution:",info$option);
  }
  return(prZ);
}#--end of function
