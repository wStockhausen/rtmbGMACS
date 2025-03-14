#'
#' @title Calculate surveys catchability for all model categories across time
#' @description
#' Function to calculate surveys catchability for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Surveys()])
#' @param params - RTMB parameters list with elements specific to survey catchability
#' @param lstSels - list of selectivity functions (output from [calcSelectivity()])
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
calcSurveys<-function(dims,info,params,lstSels,verbose=FALSE){
  if (verbose) cat("Starting calcSurvey.\n")
  lstSrvVals<-list();#--arrays of expanded survey catchability, by survey
  if (info$option=="pre-specified"){
    ##--"data" option----
    p = params$pSrv_FPs;#--vector of (fixed) selectivity values
    #--need to expand to p to all population categories, years, and seasons
    dfrUniqSrvs = info$dfrIdx2Pars |> dplyr::distinct(fcn_idx);
    for (if_ in dfrUniqSrvs$fcn_idx){
      dfrDims2Pars = info$dfrDims2Pars |> dplyr::filter(fcn_idx==if_);
      arrSrvVals = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
      for (iy_ in 1:dims$nYs){
        #--testing: iy_ = 1;
        y_ = dims$y[iy_];
        for (is_ in 1:dims$nSs){
          #--testing: is_= 1;
          s_ = dims$s[is_];
          dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |>
                       dplyr::mutate(ic_=dplyr::row_number());
          dfrIdxs = dfrDims |> dplyr::left_join(dfrDims2Pars,
                                                by = dplyr::join_by(y, s, r, x, m, p, z));
          if (nrow(dfrIdxs)>0){
            ic_ = dfrIdxs$ic_;
            arrSrvVals[iy_,is_,ic_] = p[dfrIdxs$pidx];
          }
        }#--is_ loop
      }#--iy_ loop
      lstSrvVals[[as.character(if_)]] = arrSrvVals;
      rm(arrSrvVals);
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
      } else
        stop("unrecognized selectivity function option for calcSurvey:",rwUHCs$fcn);
      }
      lstFcns[[paste(rwUHCs$fcn_idx,"+",rwUHCs$grp_idx)]] = fcn;#--save the function
      rm(fcn);
    }#--rw loop

    #--loop over surveys, years, seasons, evaluate survey catchability functions----
    for (rw in 1:nrow(info$dfrUHCs)){
      #--rw = 1;
      rwUHCs = info$dfrUHCs[rw,];
      fcn_idx = paste(rwUHCs$fcn_idx,"+",rwUHCs$grp_idx);
      fcn = lstFcns[[fcn_idx]];
      arrSrvVals = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
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
            arrSrvVals[iy_,is_,ic_] = fcn(as.numeric(dfrIdxs$z));
          }
        }#--is_ loop
      }#--iy_ loop
      lstSrvVals[[fcn_idx]] = arrSrvVals;
    }#--rw loop
  } else {
    stop("unrecognized type option for calcSurvey:",info$option);
  }
  return(lstSrvVals);
}#--end of function
