#'
#' @title Calculate fully-selected fisheries capture rates for all model categories across time
#' @description
#' Function to calculate fully-selected fisheries capture rates for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_FisheryCaptureRates()])
#' @param params - RTMB parameters list with elements specific to fisheries rates
#' @param verbose - flag to print diagnostic info
#'
#' @return A list by fishery fleet of fully-selected fishery captures as an array
#'
#' @details If `f_`, `y_`, and `s_` are the fishery fleet, year, and season of interest,
#' then the fully-selected fisheries capture rates by population category `ic_` are:
#' fshCRs = lst[[f_]][y_,s_,ic_];
#'
#' @import dplyr
#'
#' @md
#' @export
#'
calcFisheriesCaptureRates<-function(dims,info,params,verbose=FALSE){
  if (verbose) cat("Starting calcFisheriesRates.\n")
  lstFshVals<-list();#--arrays of expanded fishery capture rates, by fleet
  for (flt_ in info$flts){
    lstFshVals[[flt_]] = RTMB::AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
  }
  if (info$option=="pre-specified"){
    ##--"pre-specified" option----
    p = params$pSrvQ_FPs;#--vector of (fixed) fisheries rates values
    #--need to expand to p to all population categories, years, and seasons
    for (flt_ in info$flts){
      #--testing: flt_ = info$flts[1];
      dfrDims2Pars = info$dfrDims2Pars |> dplyr::filter(flt==flt_);
      for (iy_ in 1:dims$nYs){
        #--testing: iy_ = 1;
        y_ = dims$y[iy_];
        for (is_ in 1:dims$nSs){
          #--testing: is_= 1;
          s_ = dims$s[is_];
          dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |>
                       dplyr::mutate(ic_=dplyr::row_number());
          dfrIdxs = dfrDims |> dplyr::inner_join(dfrDims2Pars,
                                                by = dplyr::join_by(y, s, r, x, m, p, z));
          if (nrow(dfrIdxs)>0){
            ic_ = dfrIdxs$ic_;
            lstFshVals[[flt_]][iy_,is_,ic_] = p[dfrIdxs$pidx];
          }
        }#--is_ loop
      }#--iy_ loop
    }#--if_ loop
  } else if (tolower(info$option)=="function"){
    ##--"function" option----
    ###--calculate inputs to functions----
    ####--for each input parameter to a function, p = MP + OP + DP + ...
    ####--calculate values only for unique combinations of MP, OP, DP, etc.
    dfrUCs = info$dfrUniqCmbs;
    nRWs = nrow(dfrUCs);
    vals  = AD(array(0,nRWs)); #
    #browser();
    for (rw in 1:nrow(dfrUCs)){
      #--testing: rw = 1;
      dfrUCr = dfrUCs[rw,];
      p = AD(0);
      if (!is.na(dfrUCr$mpr_idx[1])) {
        p = p + params$pFCRs_MPs[dfrUCr$mpr_idx[1]];
      }
      if (any(names(dfrUCr)=="opr_idx"))
        if(!is.na(dfrUCr$opr_idx[1])) {
          if (dfrUCr$op_type=="additive") {
            p = p + params$pFCRs_OPs[dfrUCr$opr_idx[1]];
          } else {p = p * params$pFCRs_OPs[dfrUCr$opr_idx[1]];}
        }
      if (any(names(dfrUCr)=="dpr_idx"))
        if (!is.na(dfrUCr$dpr_idx[1])) {
          if (dfrUCr$dv_type=="additive") {
            p = p + params$pFCRs_DPs[dfrUCr$dpr_idx[1]];
          } else {p = p * params$pFCRs_DPs[dfrUCr$dpr_idx[1]];}
        }
      vals[rw] = p;
    }
    if (verbose){
      print(dfrUCs |> dplyr::mutate(vals_=vals));
    }

    ###--create index vector into `vals` using the index names from dfrUCs----
    idxVals = 1:length(vals);
    names(idxVals) = dfrUCs$idx;

    ###--create individual capture rate functions----
    expFCR<-function(z,pLnFCR,verbose){
      # print(pLnFCR);
      RTMB::AD(array(exp(pLnFCR),length(z)))
    }
    dZ = unname(dims$zb[2]-dims$zb[1]);#--bin size
    lstFcns = list();
    for (rw in 1:nrow(info$dfrUHCs)){
      #--for testing: rw = 1;
      rwUHCs = info$dfrUHCs[rw,];
      rwsUCs = rwUHCs |> dplyr::inner_join(dfrUCs);
      if (rwUHCs$fcn=="exp"){####--expFCR----
        fcn<-function(z){#print(vals[idxVals[rwUHCs$pLnFCR]]);
                         expFCR(z,
                                vals[idxVals[rwUHCs$pLnFCR]],
                                verbose)};
      } else {
        stop("unrecognized function option for calcFisheriesCaptureRates:",rwUHCs$fcn);
      }
      #nm = concatenateText(rwUHCs$fcn_idx,rwUHCs$grp_idx,)
      lstFcns[[rwUHCs$full_idx]] = fcn;#--save the function
      rm(fcn);
    }#--rw loop
    #browser();

    #--loop over fisheries, years, seasons, evaluate fisheries capture rates functions----
    for (rw in 1:nrow(info$dfrUHCs)){
      #--rw = 1;
      rwUHCs  = info$dfrUHCs[rw,];
      flt_    = rwUHCs$flt;
      if (verbose) cat("flt_:",flt_,"full_idx:",rwUHCs$full_idx,"\n")
      fcn = lstFcns[[rwUHCs$full_idx]];
      arrFCRsVals = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
      for (iy_ in 1:dims$nYs){
        #--iy_ = 1;
        y_ = dims$y[iy_];
        for (is_ in 1:dims$nSs){
          #--is_= 1;
          s_ = dims$s[is_];
          dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |> dplyr::mutate(ic_=dplyr::row_number());
          dfrIdxs = dfrDims |> dplyr::inner_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z)) |>
                       dplyr::filter(full_idx==rwUHCs$full_idx);
          if (nrow(dfrIdxs)>0){
            ic_ = dfrIdxs$ic_;
            arrFCRsVals[iy_,is_,ic_] = fcn(as.numeric(dfrIdxs$z));
          }
          #browser();
        }#--is_ loop
      }#--iy_ loop
      lstFshVals[[flt_]] = lstFshVals[[flt_]] + arrFCRsVals;
    }#--rw loop
  } else {
    stop("unrecognized type option for calcFisheriesCaptureRates:",info$option);
  }
  return(lstFshVals);
}#--end of function
