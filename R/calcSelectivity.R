#'
#' @title Calculate selectivity functions for all model categories across time
#' @description
#' Function to calculate the selectivity functions for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Selectivity()])
#' @param params - RTMB parameters list with elements specific to selectivity
#' @param verbose - flag to print diagnostic info
#'
#' @return A list of size-specific selectivity arrays by selectivity function
#'
#' @details If `f_`, `y_`, and `s_` are the selectivity function, year, and season of interest,
#' then the indices into the population vector and associated selectivity values are
#' given by `ic_` and `selVals`
#' selVals = lst[[f_]][y_,s_,ic_];
#'
#' @import dplyr
#'
#' @md
#' @export
#'
calcSelectivity<-function(dims,info,params,verbose=FALSE){
  if (verbose) cat("Starting calcSelectivity.\n")
  lstSelVals<-list();
  if (info$option=="data"){
    ##--"data" option----
    p = params$pSel_FPs;#--vector of selectivity values
    #--need to expand to p to all sel functions, years, and seasons
    dfrUniqSels = info$dfrIdx2Pars |> dplyr::distinct(fcn_idx);
    for (if_ in dfrUniqSels$fcn_idx){
      #--if_ = dfrUniqSels$fcn_idx[1];
      dfrDims2Pars = info$dfrDims2Pars |> dplyr::filter(fcn_idx==if_);
      arrSelVals = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
      lstYs = list();
      for (iy_ in 1:dims$nYs){
        #--iy_ = 1;
        y_ = dims$y[iy_];
        lstSs = list();
        for (is_ in 1:dims$nSs){
          #--is_= 1;
          s_ = dims$s[is_];
          dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |>
                       dplyr::mutate(ic_=dplyr::row_number());
          dfrIdxs = dfrDims |> dplyr::left_join(dfrDims2Pars,
                                                by = dplyr::join_by(y, s, r, x, m, p, z));
          if (nrow(dfrIdxs)>0){
            ic_ = dfrIdxs$ic_;
            arrSelVals[iy_,is_,ic_] = p[dfrIdxs$pidx];
            # lstSs[[names(s_)]] = list(ic_ = dfrIdxs$ic_,
            #                           sparse_idx=dfrIdxs$sparse_idx,
            #                           selVals = p[dfrIdxs$pidx]);
          }
        }#--is_ loop
        # lstYs[[names(y_)]] = lstSs;
        # rm(lstSs);
      }#--iy_ loop
      lstSelVals[[as.character(if_)]] = arrSelVals;
      # rm(lstYs);
    }#--if_ loop
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
        p = p + params$pSel_MPs[dfrUCr$mpr_idx[1]];
      }
      if (any(names(dfrUCr)=="opr_idx"))
        if(!is.na(dfrUCr$opr_idx[1])) {
          if (dfrUCr$op_type=="additive") {
            p = p + params$pSel_OPs[dfrUCr$opr_idx[1]];
          } else {p = p * params$pSel_OPs[dfrUCr$opr_idx[1]];}
        }
      if (any(names(dfrUCr)=="dpr_idx"))
        if (!is.na(dfrUCr$dpr_idx[1])) {
          if (dfrUCr$dv_type=="additive") {
            p = p + params$pSel_DPs[dfrUCr$dpr_idx[1]];
          } else {p = p * params$pSel_DPs[dfrUCr$dpr_idx[1]];}
        }
      vals[rw] = p;
    }
    if (verbose){
      print(dfrUCs |> dplyr::mutate(vals_=vals));
    }

    ###--create index vector into `vals` using the index names from dfrUCs----
    idxVals = 1:length(vals);
    names(idxVals) = dfrUCs$idx;

    ###--create individual selectivity functions----
    dZ = unname(dims$zb[2]-dims$zb[1]);#--bin size
    lstFcns = list();
    for (rw in 1:nrow(info$dfrUHCs)){
      #--for testing: rw = 1;
      rwUHCs = info$dfrUHCs[rw,];
      rwsUCs = rwUHCs |> dplyr::inner_join(dfrUCs);
      if (rwUHCs$fcn=="const_sel"){####--const_sel----
        fcn<-function(z){const_sel(z,
                                  vals[idxVals[rwUHCs$pCnst]],
                                  vals[idxVals[rwUHCs$pRefZ]],
                                  verbose=verbose)};
      } else
      if (rwUHCs$fcn=="asclogistic"){####--asclogistic----
        fcn<-function(z){asclogistic( z,
                                      vals[idxVals[rwUHCs$pZ50]],
                                      vals[idxVals[rwUHCs$pSlp]],
                                      vals[idxVals[rwUHCs$pRefZ]],
                                      verbose=verbose)};
      } else
      if (rwUHCs$fcn=="asclogistic1"){####--asclogistic1----
        fcn<-function(z){asclogistic1(z,
                                      vals[idxVals[rwUHCs$pZ50]],
                                      vals[idxVals[rwUHCs$pWdZ]],
                                      vals[idxVals[rwUHCs$pRefZ]],
                                      verbose=verbose)};
      } else
      if (rwUHCs$fcn=="asclogistic5095"){####--asclogistic5095----
        fcn<-function(z){ asclogistic5095(z,
                                          vals[idxVals[rwUHCs$pZ50]],
                                          vals[idxVals[rwUHCs$pZ95]],
                                          vals[idxVals[rwUHCs$pRefZ]],
                                          verbose=verbose)};
      } else
      if (rwUHCs$fcn=="asclogistic50D95"){####--asclogistic50D95----
        fcn<-function(z){asclogistic50D95(z,
                                          vals[idxVals[rwUHCs$pZ50]],
                                          vals[idxVals[rwUHCs$pZ9550]],
                                          vals[idxVals[rwUHCs$pRefZ]],
                                          verbose=verbose)};
      } else
      if (rwUHCs$fcn=="dbllogistic"){####--dbllogistic----
        fcn<-function(z){dbllogistic(z,
                                      vals[idxVals[rwUHCs$pAscZ50]],
                                      vals[idxVals[rwUHCs$pAscSlp]],
                                      vals[idxVals[rwUHCs$pDscZ50]],
                                      vals[idxVals[rwUHCs$pDscSlp]],
                                      vals[idxVals[rwUHCs$pRefZ]],
                                      verbose=verbose)};
      } else
      if (rwUHCs$fcn=="dbllogistic5095"){####--dbllogistic5095----
        fcn<-function(z){ dbllogistic5095(z,
                                          vals[idxVals[rwUHCs$pAscZ50]],
                                          vals[idxVals[rwUHCs$pAscZ95]],
                                          vals[idxVals[rwUHCs$pDscZ95]],
                                          vals[idxVals[rwUHCs$pDscZ50]],
                                          vals[idxVals[rwUHCs$pRefZ]],
                                          verbose=verbose)};
      } else
      if (rwUHCs$fcn=="ascnormal1"){####--ascnormal1----
        fcn<-function(z){ ascnormal1(z,
                                    vals[idxVals[rwUHCs$pAscZ1]],
                                    vals[idxVals[rwUHCs$pAscWdZ]],
                                    dZ=dZ,
                                    verbose=verbose)};
      } else
      if (rwUHCs$fcn=="ascnormal2"){####--ascnormal2----
        fcn<-function(z){ ascnormal2(z,
                                    vals[idxVals[rwUHCs$pAscZ1]],
                                    vals[idxVals[rwUHCs$pAscRefS]],
                                    vals[idxVals[rwUHCs$pRefZ]],
                                    dZ=dZ,
                                    verbose=verbose)};
      } else
      if (rwUHCs$fcn=="ascnormal2a"){####--ascnormal2a----
        fcn<-function(z){ ascnormal2a(z,
                                    vals[idxVals[rwUHCs$pAscZ1]],
                                    vals[idxVals[rwUHCs$pZatRefS]],
                                    vals[idxVals[rwUHCs$pRefs]],
                                    dZ=dZ,
                                    verbose=verbose)};
      } else
      if (rwUHCs$fcn=="ascnormal2b"){####--ascnormal2b----
        fcn<-function(z){ ascnormal2b(z,
                                    vals[idxVals[rwUHCs$pAscZ1]],
                                    vals[idxVals[rwUHCs$pDZ2RefS]],
                                    vals[idxVals[rwUHCs$pRefS]],
                                    dZ=dZ,
                                    verbose=verbose)};
      } else
      if (rwUHCs$fcn=="ascnormal3"){####--ascnormal3----
        fcn<-function(z){ ascnormal3(z,
                                    vals[idxVals[rwUHCs$pAscZ1]],
                                    vals[idxVals[rwUHCs$pSatZ2]],
                                    vals[idxVals[rwUHCs$pMxZ1]],
                                    vals[idxVals[rwUHCs$pRefZ2]],
                                    dZ=dZ,
                                    verbose=verbose)};
      } else
      if (rwUHCs$fcn=="dblnormal4"){####--dblnormal4----
        fcn<-function(z){ dblnormal4(z,
                                    vals[idxVals[rwUHCs$pAscZ1]],
                                    vals[idxVals[rwUHCs$pAscWd]],
                                    vals[idxVals[rwUHCs$pDscZ1]],
                                    vals[idxVals[rwUHCs$pDscWd]],
                                    dZ=dZ,
                                    verbose=verbose)};
      } else
      if (rwUHCs$fcn=="dblnormal4"){####--dblnormal4----
        fcn<-function(z){ dblnormal4(z,
                                    vals[idxVals[rwUHCs$pAscZ1]],
                                    vals[idxVals[rwUHCs$pAscWd]],
                                    vals[idxVals[rwUHCs$pDscZ1]],
                                    vals[idxVals[rwUHCs$pDscWd]],
                                    dZ=dZ,
                                    verbose=verbose)};
      } else
      if (rwUHCs$fcn=="dblnormal4a"){####--dblnormal4a----
        fcn<-function(z){ dblnormal4a(z,
                                    vals[idxVals[rwUHCs$pAscZ1]],
                                    vals[idxVals[rwUHCs$pAscWd]],
                                    vals[idxVals[rwUHCs$pDscSclDZ]],
                                    vals[idxVals[rwUHCs$pDscWd]],
                                    vals[idxVals[rwUHCs$pRefZ]],
                                    dZ=dZ,
                                    verbose=verbose)};
      } else
      if (rwUHCs$fcn=="dblnormal6"){####--dblnormal6----
        fcn<-function(z){ dblnormal6(z,
                                    vals[idxVals[rwUHCs$pAscZ1]],
                                    vals[idxVals[rwUHCs$pAscWd]],
                                    vals[idxVals[rwUHCs$pDscZ1]],
                                    vals[idxVals[rwUHCs$pDscWd]],
                                    vals[idxVals[rwUHCs$pAscFlr]],
                                    vals[idxVals[rwUHCs$pDscFlr]],
                                    dZ=dZ,
                                    verbose=verbose)};
      } else
      if (rwUHCs$fcn=="stackedLogistic1"){####--stackedLogistic1----
        fcn<-function(z){stackedLogistic1(z,
                                          vals[idxVals[rwUHCs$pMnZ1]],
                                          vals[idxVals[rwUHCs$pSdZ1]],
                                          vals[idxVals[rwUHCs$pMnZ2]],
                                          vals[idxVals[rwUHCs$pSdZ2]],
                                          vals[idxVals[rwUHCs$pOmga]],
                                          verbose=verbose)};
      } else {
        stop("unrecognized selectivity function option for calcSelectivity:",rwUHCs$fcn);
      }
      lstFcns[[paste(rwUHCs$fcn_idx,"+",rwUHCs$grp_idx)]] = fcn;#--save it
      rm(fcn);
    }#--rw loop

    #--loop over sel functions, years, seasons, evaluate selectivity functions----
    for (rw in 1:nrow(info$dfrUHCs)){
      #--rw = 1;
      rwUHCs = info$dfrUHCs[rw,];
      fcn_idx = paste(rwUHCs$fcn_idx,"+",rwUHCs$grp_idx);
      fcn = lstFcns[[fcn_idx]];
      arrSelVals = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
      for (iy_ in 1:dims$nYs){
        #--iy_ = 1;
        y_ = dims$y[iy_];
        for (is_ in 1:dims$nSs){
          #--is_= 1;
          s_ = dims$s[is_];
          dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |> dplyr::mutate(ic_=dplyr::row_number());
          dfrIdxs = dfrDims |> dplyr::inner_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z)) |>
                       dplyr::filter(fcn_idx==rwUHCs$fcn_idx,grp_idx==rwUHCs$grp_idx);
          if (nrow(dfrIdxs)>0){
            ic_ = dfrIdxs$ic_;
            arrSelVals[iy_,is_,ic_] = fcn(as.numeric(dfrIdxs$z));
          }
        }#--is_ loop
      }#--iy_ loop
      lstSelVals[[fcn_idx]] = arrSelVals;
    }#--rw loop
  } else {
    stop("unrecognized type option for calcSelectivity:",info$option);
  }
  return(lstSelVals);
}#--end of function
