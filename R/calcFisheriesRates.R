#'
#' @title Calculate size-specific fisheries capture and mortality rates for all model categories across time
#' @description
#' Function to calculate sie-specific fisheries capture and mortality rates for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_FisheriesRates()])
#' @param params - RTMB parameters list with elements specific to survey catchability
#' @param lstFCRs - list of selectivity value arrays (output from [calcFisheriesCaptureRates()])
#' @param lstHMRs - list of selectivity value arrays (output from [calcFisheriesHandlingMortalityRates()])
#' @param lstSels - list of selectivity value arrays (output from [calcSelectivity()])
#' @param verbose - flag to print diagnostic info
#'
#' @return A list by fishery fleet of size-specific fishery capture, retention, discard,
#' and total mortality rates as an array
#'
#' @details If `f_`, `y_`, and `s_` are the fishery fleet, year, and season of interest,
#' then the indices into the population vector and associated fishing rates are
#' given by `ic_` and `type_`:
#'
#' fshRs = lst[[f_]][type_,y_,s_,ic_];
#'
#' where `type_` is one of "FCR" (fishery capture rate),"RMR" (retention mortality rate),
#' "DMR" (discard mortality rate), "TFMR" (total fishing mortality rate)
#'
#' @import dplyr
#'
#' @md
#' @export
#'
calcFisheriesRates<-function(dims,info,lstFCRs,lstHMRs,lstSels,verbose=FALSE){
  if (verbose) cat("Starting calcFisheriesRates.\n")
  nTypes = 8; #--"FCR","RMR","DMR","TFMR"
  lstFshVals<-list();#--arrays of expanded fishery rates, by fleet
  for (flt_ in unique(info$dfrUniqCmbs$flt)){
    lstFshVals[[flt_]] = RTMB::AD(array(0,c(nTypes,dims$nYs,dims$nSs,dims$nCs)));
  }

  #--loop over unique combinations of fleet, capture selectivity, retention curve----
  for (rw in 1:nrow(info$dfrUniqCmbs)){
    #--rw = 1;
    rwUCs  = info$dfrUniqCmbs[rw,];
    flt_    = rwUCs$flt;
    cp_idx = rwUCs$cap_idx;#--capture selectivity index
    rt_idx = rwUCs$ret_idx;#--retention curve index
    if (verbose) cat("flt_:",flt_,"cap_idx:",cp_idx,"ret_idx:",rt_idx,"\n")
    arrCapVals = lstSels[[cp_idx]];
    if (rwUCs$ret_idx>0) {#--retention curve
      arrRetVals = lstSels[[rt_idx]];
    } else {
      arrRetVals = RTMB::AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
    }
    arrFCRs = lstFCRs[[flt_]]*arrCapVals;                     #--capture rates
    arrRMRs = arrFCRs * arrRetVals;                           #--retained mortality rates
    arrDMRs = arrFCRs * lstHMRs[[flt_]] * (AD(1)-arrRetVals); #--discard mortality rates
    for (iy_ in 1:dims$nYs){
      #--iy_ = 1;
      y_ = dims$y[iy_];
      for (is_ in 1:dims$nSs){
        #--is_= 1;
        s_ = dims$s[is_];
        dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |> dplyr::mutate(ic_=dplyr::row_number());
        dfrIdxs = dfrDims |>
                    dplyr::inner_join(info$dfrCmbs |> dplyr::filter(flt==flt_,cap_idx==cp_idx,ret_idx==rt_idx),
                                      by = dplyr::join_by(y, s, r, x, m, p, z));
        if (nrow(dfrIdxs)>0){
          ic_ = dfrIdxs$ic_;
          lstFshVals[[flt_]][1,iy_,is_,ic_] = lstFshVals[[flt_]][1,iy_,is_,ic_] + arrFCRs[iy_,is_,ic_];                     #--capture rates
          lstFshVals[[flt_]][2,iy_,is_,ic_] = lstFshVals[[flt_]][2,iy_,is_,ic_] + arrRMRs[iy_,is_,ic_];                     #--retained mortality rates
          lstFshVals[[flt_]][3,iy_,is_,ic_] = lstFshVals[[flt_]][3,iy_,is_,ic_] + arrDMRs[iy_,is_,ic_];                     #--discard mortality rates
          lstFshVals[[flt_]][4,iy_,is_,ic_] = lstFshVals[[flt_]][4,iy_,is_,ic_] + arrRMRs[iy_,is_,ic_]+arrDMRs[iy_,is_,ic_];#--total fishing mortality rates
          lstFshVals[[flt_]][5,iy_,is_,ic_] = lstFshVals[[flt_]][5,iy_,is_,ic_] + lstFCRs[[flt_]][iy_,is_,ic_];             #--fully-selected capture rates
          lstFshVals[[flt_]][6,iy_,is_,ic_] = lstFshVals[[flt_]][6,iy_,is_,ic_] + lstHMRs[[flt_]][iy_,is_,ic_];             #--handling mortality rates
          lstFshVals[[flt_]][7,iy_,is_,ic_] = lstFshVals[[flt_]][7,iy_,is_,ic_] + arrCapVals[iy_,is_,ic_];                  #--capture selectivity
          lstFshVals[[flt_]][8,iy_,is_,ic_] = lstFshVals[[flt_]][8,iy_,is_,ic_] + arrRetVals[iy_,is_,ic_];                  #--retention function
        }
      }#--is_ loop
    }#--iy_ loop
  }#--rw loop
  return(lstFshVals);
}#--end of function
