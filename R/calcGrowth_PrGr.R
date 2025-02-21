#'
#' @title Calculate growth
#' @description Function to calculate growth.
#' @param pGrA - mean post-molt size at zSclGrA
#' @param zGrA - pre-molt size at yielding pGrA as mean post-molt size
#' @param pGrB - mean post-molt size at zSclGrB
#' @param zGrB - pre-molt size at yielding pGrB as mean post-molt size
#' @param pGrBeta - gamma distribution scale parameter for post-molt variability
#' @param zBs_from - pre-molt sizes at which to calculate values
#' @param zBs_to - post-molt sizes at which to calculate values
#' @return growth matrix
#'
#' @details The formula used is
#'
#' mnZs = pGrA\*exp(log(pGrB/pGrA)/log(zGrB/zGrA)\*log(zBs/zGrA));
#'
#' @examples
#' # example code
#' grM = grwPwrLaw1(33,25,150,125,1,seq(25,180,5),seq(25,180,5));
#'
#' @md
#' @export
#'
grwPwrLaw1<-function(pGrA,zGrA,pGrB,zGrB,pGrBeta,zBs_from,zBs_to){
  mnZs = pGrA*exp(log(pGrB/pGrA)/log(zGrB/zGrA)*log(zBs_from/zGrA));
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
grwPwrLaw2<-function(pGrA,pGrB,pGrBeta,zBs){
  mnZs = exp(grA+grB*log(zBs));
  return(grM);
}

#'
#' @title Calculate the probability of undergoing transitions between size categories during growth, for all model categories across time
#' @description
#' Function to calculate the probability of undergoing transitions between size categories during growth, for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Growth_PrGr()])
#' @param params - RTMB parameters list with elements specific to the probability of undergoing a transition between size categories
#' @param verbose - flag to print diagnostic info
#'
#' @return TODO: might want to return a list of a list of matrices
#'
#' @details TBD
#'
#' @import dplyr
#'
#' @md
#' @export
#'
calcGrowth_PrGr<-function(dims,info,params,verbose=FALSE){
  if (verbose) cat("Starting calcGrowth_PrGr.\n")
  diags = array(c(1:dims$nCs,1:dims$nCs),dim=c(dims$nCs,2));
  lstY = list();
  if (info$option=="data"){
    ##--"data" option----
    p_ = params$pPrGr_FPs;#--vector of weights-at-size
    #--need to expand to p to all years, seasons, and population categories
    for (iy_ in 1:dims$nYs){
      #--iy_ = 2;
      y_ = dims$y[iy_];
      if (verbose) cat("year: ",y_,"\n")
      lstS = list();
      for (is_ in 1:dims$nSs){
        #--is_= 1;
        s_ = dims$s[is_];
        if (verbose) cat("season: ",s_,"\n")
        tmPrGr = AD(array(0,c(dims$nCs,dims$nCs)));#--rows: `to` category; columns: `from` category
        dfrDims = dims$dmsYSC |> dplyr::filter(y==y_,s==s_);
        idx_base = dfrDims$sparse_idx[1]-1;
        dfrDimsTF = (dfrDims |> dplyr::rename(z_from=z,col_idx=sparse_idx)) |>
                      dplyr::left_join((dfrDims |> dplyr::rename(z_to=z,row_idx=sparse_idx)),
                                       by=dplyr::join_by(y, s, r, x, m, p),
                                       relationship="many-to-many") |>
                     dplyr::select(y,s,r,x,m,p,z_from,z_to,col_idx,row_idx) |>
                     dplyr::mutate(col_idx=col_idx-idx_base,
                                   row_idx=row_idx-idx_base);
        #--expand
        dfrIdxs = dfrDimsTF |>
                     dplyr::left_join(info$dfrDims2Pars,
                                      by =dplyr::join_by(y, s, r, x, m, p, z_from,z_to)) |>
                     dplyr::filter(!is.na(pidx));
        #--assign defined transition probabilities
        tmPrGr[as.matrix(dfrIdxs[,c("row_idx","col_idx")])] = p_[dfrIdxs$pidx];
        tmPrGr[diags] = AD(0);                 #--remove defined self-transition probabilities
        tmPrGr[diags] = AD(1)-colSums(tmPrGr);#--assign self-transition probabilities
        lstS[[names(s_)]] = tmPrGr;
      }#--is_ loop
      lstY[[names(y_)]] = lstS;
    }#--iy_ loop
  } else if (tolower(info$option)=="function"){
    ##--"function" option----
    ###--calculate inputs to functions----
    ####--for each input parameter to a function, p = MP + OP + DP + ...
    ####--calculate values only for unique combinations of MP, OP, DP, etc.
    if (verbose) {
      cat("function option.\n")
      cat("\nCalculating inputs to functions.\n")
      print(params$pPrGr_MPs);
    }

    dfrUCs = info$dfrUniqCmbs;
    nRWs = nrow(dfrUCs);
    vals  = AD(array(0,nRWs)); #
    for (rw in 1:nrow(dfrUCs)){
      #--testing: rw = 1;
      dfrUCr = dfrUCs[rw,];
      # print(dfrUCr);
      p = AD(0);
      if (!is.na(dfrUCr$mpr_idx[1])) {
        p = p + params$pPrGr_MPs[dfrUCr$mpr_idx[1]];
        # print(p);
      }
      if (!is.null(dfrUCr$opr_idx))
        if (!is.na(dfrUCr$opr_idx[1])){
          if (dfrUCr$op_type=="additive") {
            p = p + params$pPrGr_OPs[dfrUCr$opr_idx[1]];
          } else {p = p * params$pPrGr_OPs[dfrUCr$opr_idx[1]];}
        }
      if (!is.null(dfrUCr$dpr_idx))
        if (!is.na(dfrUCr$dpr_idx[1])) {
          if (dfrUCr$dv_type=="additive") {
            p = p + params$pPrGr_DPs[dfrUCr$dpr_idx[1]];
          } else {p = p * params$pPrGr_DPs[dfrUCr$dpr_idx[1]];}
        }
      #--TBD: add in missing components (REs, env covars, etc)
      # print(p);
      vals[rw] = p;
    }
    if (verbose){
      print(vals);
      #print(dfrUCs |> dplyr::mutate(vals_=as.vector(vals)));
    }

    ###--create index vector into `vals` using the index names from dfrUCs----
    idxVals = 1:length(vals);
    names(idxVals) = dfrUCs$idx;

    ###--calculate transition matrices for growth----
    if (verbose) cat("Calculating transition matrices for growth.\n")
    lstY = list();
    for (iy_ in 1:dims$nYs){
      #--iy_ = 1;
      y_ = dims$y[iy_];
      if (verbose) cat("year: ",y_,"\n")
      lstS = list();
      for (is_ in 1:dims$nSs){
        #--is_= 1;
        s_ = dims$s[is_];
        if (verbose) cat("season: ",s_,"\n")
        tmPrGr = AD(array(0,c(dims$nCs,dims$nCs)));#--rows: `to` category; columns: `from` category
        dfrDims  = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |> dplyr::mutate(ic_=dplyr::row_number());
        idx_base = dfrDims$sparse_idx[1]-1;
        dfrDimsTF = (dfrDims |> dplyr::rename(z_from=z,col_idx=sparse_idx)) |>
                      dplyr::left_join((dfrDims |> dplyr::rename(z_to=z,row_idx=sparse_idx)),
                                       by=dplyr::join_by(y, s, r, x, m, p),
                                       relationship="many-to-many") |>
                     dplyr::select(y,s,r,x,m,p,z_from,z_to,col_idx,row_idx) |>
                     dplyr::mutate(col_idx=col_idx-idx_base,
                                   row_idx=row_idx-idx_base);

        dfrIdxsA = dfrDimsTF |> dplyr::left_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z_from, z_to));
        dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)==tolower("grwPwrLaw1"),as.numeric(z_from)<=as.numeric(z_to));
        if ((nRWs=nrow(dfrIdxs)) > 0){
          #--parameters are: pGrA,zGrA,pGrB,zGrB,pGrBeta
          col_idxs = dfrIdxs$col_idx;#--from
          row_idxs = dfrIdxs$row_idx;#--to
          tmPrGr[array(c(row_idxs,col_idxs),dim=c(nRWs,2))] = grwPwrLaw1(vals[idxVals[dfrIdxs$pGrA]],
                                                                         vals[idxVals[dfrIdxs$zGrA]],
                                                                         vals[idxVals[dfrIdxs$pGrB]],
                                                                         vals[idxVals[dfrIdxs$zGrB]],
                                                                         vals[idxVals[dfrIdxs$pGrBeta]],
                                                                         dfrIdxs$z_from,dfrIdxs$z_to);
        }
        dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)=="grwPwrLaw2");
        if ((nRWs=nrow(dfrIdxs)) > 0){
          #--parameters are: ????
          col_idxs = dfrIdxs$col_idx;#--from
          row_idxs = dfrIdxs$row_idx;#--to
          tmPrGr[array(c(row_idxs,col_idxs),dim=c(nRWs,2))] = grwPwrLaw2(vals[idxVals[dfrIdxs$pGrA]],
                                                                         vals[idxVals[dfrIdxs$pGrB]],
                                                                         vals[idxVals[dfrIdxs$pGrBeta]],
                                                                         dfrIdxs$z);
        }
        tmPrGr[diags] = AD(0);                 #--remove defined self-transition probabilities
        tmPrGr[diags] = AD(1)-colSums(tmPrGr);#--assign self-transition probabilities
        lstS[[names(s_)]] = tmPrGr;
      }#--is_ loop
      lstY[[names(y_)]] = lstS;
    }#--iy_ loop
  } else {
    stop("unrecognized type option for calcGrowth_PrGr:",info$option);
  }
  if (verbose) cat("Finished calcGrowth_PrGr().\n")
  return(lstY);
}#--end of function
