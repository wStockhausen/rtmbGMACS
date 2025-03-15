#'
#' @title Calculate surveys catchability for all model categories across time
#' @description
#' Function to calculate surveys catchability for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_SurveysCatchability()])
#' @param params - RTMB parameters list with elements specific to survey catchability
#' @param lstSels - list of selectivity value arrays (output from [calcSelectivity()])
#' @param verbose - flag to print diagnostic info
#'
#' @return A list of size-specific survey catchability arrays by survey fleet
#'
#' @details If `f_`, `y_`, and `s_` are the survey fleet, year, and season of interest,
#' then the indices into the population vector and associated selectivity values are
#' given by `ic_` and `selVals`
#' srvQs = lst[[f_]][y_,s_,ic_];
#'
#' @import dplyr
#'
#' @md
#' @export
#'
calcSurveysCatchability<-function(dims,info,params,lstSels,verbose=FALSE){
  if (verbose) cat("Starting calcSurveyCatchability.\n")
  lstSrvVals<-list();#--arrays of expanded survey catchability, by survey
  for (flt_ in info$flts){
    lstSrvVals[[flt_]] = RTMB::AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
  }
  if (info$option=="pre-specified"){
    ##--"data" option----
    p = params$pSrvQ_FPs;#--vector of (fixed) selectivity values
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
            lstSrvVals[[flt_]][iy_,is_,ic_] = p[dfrIdxs$pidx];
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
    for (rw in 1:nrow(dfrUCs)){
      #--testing: rw = 1;
      dfrUCr = dfrUCs[rw,];
      p = AD(0);
      if (!is.na(dfrUCr$mpr_idx[1])) {
        p = p + params$pSrv_MPs[dfrUCr$mpr_idx[1]];
      }
      if (any(names(dfrUCr)=="opr_idx"))
        if(!is.na(dfrUCr$opr_idx[1])) {
          if (dfrUCr$op_type=="additive") {
            p = p + params$pSrv_OPs[dfrUCr$opr_idx[1]];
          } else {p = p * params$pSrv_OPs[dfrUCr$opr_idx[1]];}
        }
      if (any(names(dfrUCr)=="dpr_idx"))
        if (!is.na(dfrUCr$dpr_idx[1])) {
          if (dfrUCr$dv_type=="additive") {
            p = p + params$pSrv_DPs[dfrUCr$dpr_idx[1]];
          } else {p = p * params$pSrv_DPs[dfrUCr$dpr_idx[1]];}
        }
      vals[rw] = p;
    }
    if (verbose){
      print(dfrUCs |> dplyr::mutate(vals_=vals));
    }

    ###--create index vector into `vals` using the index names from dfrUCs----
    idxVals = 1:length(vals);
    names(idxVals) = dfrUCs$idx;

    ###--create individual catchability functions----
    dZ = unname(dims$zb[2]-dims$zb[1]);#--bin size
    lstFcns = list();
    for (rw in 1:nrow(info$dfrUHCs)){
      #--for testing: rw = 1;
      rwUHCs = info$dfrUHCs[rw,];
      rwsUCs = rwUHCs |> dplyr::inner_join(dfrUCs);
      if (rwUHCs$fcn=="lnQ"){####--lnQ----
        fcn<-function(z){lnQ(z,
                             vals[idxVals[rwUHCs$pLnQ]],
                             verbose=verbose)};
      } else {
        stop("unrecognized selectivity function option for calcSurveysCatchability:",rwUHCs$fcn);
      }
      lstFcns[[paste(rwUHCs$fcn_idx,"+",rwUHCs$grp_idx)]] = fcn;#--save the function
      rm(fcn);
    }#--rw loop
    #browser();

    #--loop over surveys, years, seasons, evaluate survey catchability functions----
    for (rw in 1:nrow(info$dfrUHCs)){
      #--rw = 1;
      rwUHCs  = info$dfrUHCs[rw,];
      flt_    = rwUHCs$flt;
      fcn_idx = paste(rwUHCs$fcn_idx,"+",rwUHCs$grp_idx);
      if (verbose) cat("flt_:",flt_,"fcn_idx:",fcn_idx,"sel_idx:",rwUHCs$sel_idx,"avl_idx:",rwUHCs$avl_idx,"\n")
      arrSelVals = lstSels[[rwUHCs$sel_idx]];
      if (rwUHCs$avl_idx>0) {#--apply availability
        arrSelVals = lstSels[[rwUHCs$avl_idx]] * arrSelVals;
      }
      fcn = lstFcns[[fcn_idx]];
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
            if (rwUHCs$sel_idx>0) {
              #--combine fully-selected Q and selectivity/availability values
              lstSrvVals[[flt_]][iy_,is_,ic_] = fcn(as.numeric(dfrIdxs$z))*arrSelVals[iy_,is_,ic_];
            } else {
              #--
              lstSrvVals[[flt_]][iy_,is_,ic_] = fcn(as.numeric(dfrIdxs$z));
            }
          }
        }#--is_ loop
      }#--iy_ loop
    }#--rw loop
  } else {
    stop("unrecognized type option for calcSurveysCatchability:",info$option);
  }
  return(lstSrvVals);
}#--end of function
