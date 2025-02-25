#'
#' @title Calculate a gamma-distrbuted size istribution at recruitment
#' @description Function to calculate a gamma-distrbuted size istribution at recruitment.
#' @param mnZ - mean size at recruitment
#' @param wdZ - width (std. dev of gamma distribution)
#' @param zMn - minimum size with non-zero probability of recruitment
#' @param zMx - maximum size with non-zero probability of recruitment
#' @param zBs - sizes at which to calculate vector
#' @param dZ - size bin width to use to set
#' @return object with same dimensions as `zBs`.
#'
#' @details
#'
#' The formula used is
#'
#' $$recTS(zBs < zMn or zBs > zMx)   = 0.0$$
#' otherwise, for mnZ <=zBs<=zMx,
#' $$recTSmn = pgamma(mnZ-dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ)$$
#' $$recTSmx = pgamma(zMx+dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ)$$
#' $$recTS(zBs) = pgamma(zBs+dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ) -
#'              pgamma(zBs-dZ/2,shape=(mnZ^2)/(wdZ^2),scale=(wdZ^2)/mnZ)$$
#' $$recTS(zBs) = recTS(zBs \le zMx)/(recTSMx-recTSmn)$$
#'
#' @examples
#' # example code
#' zBs = seq(25,100,5);
#' dZ  = z[2]-z[1];
#' recTS = prRecZ1(35,20,55,25,80,zBs,dZ);
#'
#' @md
#' @export
#'
# prRecZ1<-function(mnZ,wdZ,zMn,zMx,zBs,dZ){
#   recTS = AD(array(0,dim=length(zBs)));
#   zBs = as.numeric(zBs);
#   print(zBs);
#   print(mnZ);
#   print(wdZ);
#   print(zMn);
#   print(zMx);
#   dzMn = as.double(zMn);#--need to convert to double for comparison (doesn't seem to work)
#   print(dzMn);
#   dzMx = as.double(zMx);#--need to convert to double for comparison (doesn't seem to work)
#   print(dzMx);
#   iZBs = which((dzMn<=zBs)&(zBs <= dzMx));
#   print(iZBs);
#   print(zBs[iZBs])
#   shp = (mnZ[iZBs]^2)/(wdZ[iZBs]^2);
#   print(shp);
#   scl = (wdZ[iZBs]^2)/mnZ[iZBs];
#   print(scl);
#   recTSmn = pgamma(zMn[iZBs]-dZ/2,shape=shp,scale=scl);
#   recTSmx = pgamma(zMx[iZBs]+dZ/2,shape=shp,scale=scl);
#   recTS[iZBs] = pgamma(zBs[iZBs]+dZ/2,shape=shp,scale=scl) -
#               pgamma(zBs[iZBs]-dZ/2,shape=shp,scale=scl);
#   print(recTS[iZBs]);
#   recTS[iZBs] = recTS[iZBs]/(recTSmx- recTSmn);
#   print(recTS[iZBs]);
#   print(sum(recTS[iZBs]));
#   # cat("\t",recTS,"\n");
#   return(recTS);
# }
prRecZ1<-function(mnZ,wdZ,zMn,zMx,zBs,dZ){
  recTS = AD(array(0,dim=length(zBs)));
  zBs = as.numeric(zBs);
  # cat("dZ:",dZ,"\n");
  # cat("ZBs\n");
  # print(zBs);
  # cat("mnZ\n");
  # print(mnZ);
  # cat("wdZ\n");
  # print(wdZ);
  # cat("zMn\n");
  # print(zMn);
  # cat("zMx\n");
  # print(zMx);
  #--calculate squarewave to "box" size range
  sqw = squarewave(zMn,zMx,zBs,dZ)
  # cat("sqw\n");
  # print(sqw);
  # cat("1-sqw\n");
  # print(1-sqw);
  shp = (mnZ^2)/(wdZ^2);
  # cat("shp\n");
  # print(shp);
  scl = (wdZ^2)/mnZ;
  # cat("scl\n");
  # print(scl);
  recTSmn = pgamma(zMn-dZ/2,shape=shp,scale=scl);
  recTSmx = pgamma(zMx+dZ/2,shape=shp,scale=scl);
  recTS   = pgamma(zBs+dZ/2,shape=shp,scale=scl) -
            pgamma(zBs-dZ/2,shape=shp,scale=scl);
  # cat("recTS before scaling\n")
  # print(recTS);
  # print(recTSmx- recTSmn);
  #--normalize and apply squarewave window so sum  across zBs is 1 for each population category
  recTS = sqw*recTS/(recTSmx- recTSmn);
  # cat("recTS after scaling\n")
  # print(recTS);
  # cat("sum(recTS):",sum(recTS),"\n");
  return(recTS);
}


#'
#' @title Calculate the bulk recruitment time series across time
#' @description
#' Function to calculate the bulk recruitment time series across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Recruitment_TimeSeries()])
#' @param params - RTMB parameters list with elements specific to the recruitment time series
#' @param verbose - flag to print diagnostic info
#'
#' @return TODO: might want to return a list of a list of matrices
#'
#' @details Combination of the recruitment time series, with other aspects of
#' recruitment for population category proportions and sizes, yields
#' recruitment by year, season, population category (excluding size) and size.
#'
#' @import dplyr
#'
#' @md
#' @export
#'
calcRecruitment_TimeSeries<-function(dims,info,params,verbose=FALSE,loopIC_=FALSE){
  if (verbose) cat("Starting calcRecruitment_TimeSeries.\n")
  recTS = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
  if (info$option=="data"){
    ##--"data" option----
    p = params$pRecTS_FPs;#--vector of recruitment time series values
    #--need to expand to p to all years, seasons, and population categories
    for (iy_ in 1:dims$nYs){
      #--iy_ = 1;
      y_ = dims$y[iy_];
      for (is_ in 1:dims$nSs){
        #--is_= 1;
        s_ = dims$s[is_];
        dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_));
        dfrIdxs = dfrDims |> dplyr::left_join(info$dfrDims2Pars,
                                              by = dplyr::join_by(y, s, r, x, m, p, z));
        pidx = dfrIdxs$pidx;
        recTS[iy_,is_,] = p[pidx];
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
                recTS[iy_,is_,ic_] = prRecZ1(vals[idxVals[dfrIdxs$pMnZ]],
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
            recTS[iy_,is_,ic_] = prRecZ1(vals[idxVals[dfrIdxs$pMnZ]],
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
    stop("unrecognized type option for calcRecruitment_TimeSeries:",info$option);
  }
  return(recTS);
}#--end of function
