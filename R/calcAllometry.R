#'
#' @title Calculate allometry for all model categories across time
#' @description
#' Function to calculate allometryfor all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Allometry()])
#' @param params - RTMB parameters list with allometry-specific elements
#' @param verbose - flag to print diagnostic info
#'
#' @return list
#'
#' @export
#'
calcAllometry<-function(dims,info,params,verbose=FALSE){
  if (verbose) cat("Starting calcAllometry.\n")
  wAtZ = AD(array(0,c(dims$nYs,dims$nSs,dims$nCs)));
  if (info$option=="data"){
    ##--"data" option----
    p = params$pAllom_FPs;#--vector of weights-at-size
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
        #   wAtZ[iy_,is_,ic_] = p[pidx];
        # }
        dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_));
        dfrIdxs = dfrDims |> dplyr::left_join(info$dfrDims2Pars,
                                              by = join_by(y, s, r, x, m, p, z));
        pidx = dfrIdxs$pidx;
        wAtZ[iy_,is_,] = p[pidx];
      }#--is_ loop
    }#--iy_ loop
  } else if (info$option=="function"){
    ##--"function" option----
    dfrFcns = info$Fcns$dfrIdxs |> dplyr::select(fnc,
                                                 fcn_idx);
    dfrMPs  = info$MPs$dfrMP3s |> dplyr::select(param,
                                                fcn_idx,
                                                mp_inp_idx=inp_idx,
                                                mp_pv_idx=pv_idx);
    dfrOPs  = info$OPs$dfrOP3s |> dplyr::select(param,
                                                mp_inp_idx=inp_par_idx,
                                                op_inp_idx=inp_off_idx,
                                                op_pv_idx=pv_idx,
                                                op_type=offset_type);
    dfrDPs  = info$DPs$dfrDP4s |> dplyr::select(param,
                                                mp_inp_idx=inp_par_idx,
                                                dp_inp_idx=inp_dvec_idx,
                                                dp_pv_idx=dev_par_idx,
                                                dp_dev_type=dev_type);
    dfrCmb = dfrFcns |>
               dplyr::left_join(dfrMPs) |>
               dplyr::left_join(dfrOPs) |>
               dplyr::left_join(dfrDPs,relationship="many-to-many");
    for (iy_ in 1:dims$nYs){
      #--iy_ = 1;
      y_ = dims$y[iy_];
      for (is_ in 1:dims$nSs){
        #--is_= 1;
        s_ = dims$s[is_];
        for (ic_ in 1:dims$nCs){
          #--ic_ = 1;
          dfrDims = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_))[ic_,];
          dfrIdxs = dfrDims |> dplyr::left_join(info$dfrDims2Pars);
          pidx = dfrIdxs$pidx[1];
          wAtZ[iy_,is_,ic_] = p[pidx];
        }
      }#--is_ loop
    }#--iy_ loop
  } else {
    stop("unrecognized type option for allometry:",info$option);
  }
  return(wAtZ);
}#--end of function
