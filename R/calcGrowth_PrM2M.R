#'
#' @title Calculate a size-dependent, ascending logistic function for the probability of molting to maturity
#' @description Function to calculate a size-dependent growth, ascending logistic function for the probability of molting to maturity.
#' @param mdZ - size at which curve starts to descend from 1
#' @param wdZ - width of descent (standard deviation of normal curve)
#' @param zBs - pre-molt sizes at which to calculate vector
#' @param dZ - size bin width to use to set
#' @return object with same dimensions as `zBs`.
#'
#' @details
#'
#' The formula used is
#'
#' $$prM(zBs) = exp(zBs-mdZ) \over (1+exp(zBs-mdZ))$$
#'
#' @examples
#' # example code
#' z = seq(25,100,5);
#' prM2M = prM2M_AscLogistic(55,30,z);
#'
#' @md
#' @export
#'
prM2M_AscLogistic<-function(mdZ,wdZ,zBs){
  # cat("In AscLogistic\n")
  # print(mdZ);
  # print(wdZ);
  # print(zBs);
  prM2M = AD(array(0,dim=length(zBs)));
  zBs = as.numeric(zBs);
  prM2M = AD(exp(zBs-mdZ)/(1+exp(zBs-mdZ)));
  return(prM2M);
}

#'
#' @title Calculate a size-dependent, ascending normal probability of molting to maturity
#' @description Function to calculate a size-dependent growth, descending normal probability of molting to maturity.
#' @param mdZ - size at which curve ascends to 1
#' @param wdZ - width of ascent (standard deviation of normal curve)
#' @param zBs - pre-molt sizes at which to calculate vector
#' @param dZ - size bin width to use to set
#' @return object with same dimensions as `zBs`.
#'
#' @details
#'
#' The formula used is
#'
#' $$prM(zBs > mdZ) = 1.0$$
#' $$prM(zBs \le mdZ) = exp(0.5*((zBs-mdZ) \over wdZ)^2)$$
#'
#' @examples
#' # example code
#' z = seq(25,100,5);
#' prM2M = prM2M_AscNormal(55,30,z);
#'
#' @md
#' @export
#'
prM2M_AscNormal<-function(mdZ,wdZ,zBs){
  prM2M = AD(array(0,dim=length(zBs)));
  zBs = as.numeric(zBs);
  # cat("in prM2M_AscNormal\n")
  # cat("\t",zBs," \n");
  # cat("\t",exp(-0.5*((zBs-mdZ)/wdZ)^2)," \n");
  # cat("\t",squarewave_left(mdZ,zBs),"\n")
  prM2M = AD(exp(-0.5*((zBs-mdZ)/wdZ)^2)*squarewave_left(mdZ,zBs)) +
           squarewave_right(mdZ,zBs);
  # cat("\t",prM,"\n");
  return(prM2M);
}

#'
#' @title Calculate the probability of undergoing a transition between maturity categories for all model categories across time
#' @description
#' Function to calculate the probability of undergoing a transition between maturity categories for all model categories across time.
#' @param dims - dimensions list
#' @param info - info list (output list from [extractParamInfo_Growth_PrM2M()])
#' @param params - RTMB parameters list with elements specific to the probability of undergoing a transition between maturity categories
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
calcGrowth_PrM2M<-function(dims,info,params,verbose=FALSE){
  if (verbose) cat("Starting calcGrowth_PrM2M.\n")
  diags = array(c(1:dims$nCs,1:dims$nCs),dim=c(dims$nCs,2));
  lstY = list();
  if (info$option=="data"){
    ##--"data" option----
    p_ = params$pPrM2M_FPs;#--vector of weights-at-size
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
        tmPrM2M = AD(array(0,c(dims$nCs,dims$nCs)));#--rows: `to` category; columns: `from` category
        dfrDims = dims$dmsYSC |> dplyr::filter(y==y_,s==s_);
        idx_base = dfrDims$sparse_idx[1]-1;
        dfrDimsTF = (dfrDims |> dplyr::rename(m_from=m,col_idx=sparse_idx)) |>
                      dplyr::left_join((dfrDims |> dplyr::rename(m_to=m,row_idx=sparse_idx)),
                                       by=dplyr::join_by(y, s, r, x, p, z),
                                       relationship="many-to-many") |>
                     dplyr::select(y,s,r,x,p,z,m_from,m_to,col_idx,row_idx) |>
                     dplyr::mutate(col_idx=col_idx-idx_base,
                                   row_idx=row_idx-idx_base);
        #--expand
        dfrIdxs = dfrDimsTF |>
                     dplyr::left_join(info$dfrDims2Pars,
                                      by =dplyr::join_by(y, s, r, x, p, z,m_from,m_to)) |>
                     dplyr::filter(!is.na(pidx));
        #--assign defined transition probabilities
        tmPrM2M[as.matrix(dfrIdxs[,c("row_idx","col_idx")])] = p_[dfrIdxs$pidx];
        tmPrM2M[diags] = AD(0);                 #--remove defined self-transition probabilities
        tmPrM2M[diags] = AD(1)-colSums(tmPrM2M);#--assign self-transition probabilities
        lstS[[names(s_)]] = tmPrM2M;
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
      print(params$pPrM2M_MPs)
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
        p = p + params$pPrM2M_MPs[dfrUCr$mpr_idx[1]];
        # print(p);
      }
      if (!is.null(dfrUCr$opr_idx))
        if (!is.na(dfrUCr$opr_idx[1])){
          if (dfrUCr$op_type=="additive") {
            p = p + params$pPrM2M_OPs[dfrUCr$opr_idx[1]];
          } else {p = p * params$pPrM2M_OPs[dfrUCr$opr_idx[1]];}
        }
      if (!is.null(dfrUCr$dpr_idx))
        if (!is.na(dfrUCr$dpr_idx[1])) {
          if (dfrUCr$dv_type=="additive") {
            p = p + params$pPrM2M_DPs[dfrUCr$dpr_idx[1]];
          } else {p = p * params$pPrM2M_DPs[dfrUCr$dpr_idx[1]];}
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

    ###--calculate transition matrices for molt to maturity----
    if (verbose) cat("Calculating transition matrices molt to maturity.\n")
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
        tmPrM2M = AD(array(0,c(dims$nCs,dims$nCs)));#--rows: `to` category; columns: `from` category
        dfrDims  = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |> dplyr::mutate(ic_=dplyr::row_number());
        idx_base = dfrDims$sparse_idx[1]-1;
        dfrDimsTF = (dfrDims |> dplyr::rename(m_from=m,col_idx=sparse_idx)) |>
                      dplyr::left_join((dfrDims |> dplyr::rename(m_to=m,row_idx=sparse_idx)),
                                       by=dplyr::join_by(y, s, r, x, p, z),
                                       relationship="many-to-many") |>
                     dplyr::select(y,s,r,x,p,z,m_from,m_to,col_idx,row_idx) |>
                     dplyr::mutate(col_idx=col_idx-idx_base,
                                   row_idx=row_idx-idx_base);

        dfrIdxsA = dfrDimsTF |> dplyr::left_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, p, z, m_from, m_to));
        dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)=="asclogistic");
        if ((nRWs=nrow(dfrIdxs)) > 0){
          #--inputs are pZ50, pWdZ, z
          col_idxs = dfrIdxs$col_idx;#--from
          row_idxs = dfrIdxs$row_idx;#--to
          # cat("Printing dfrIdxs$pMdZ indices to pZ50\n");
          # print(dfrIdxs$pZ50)
          # cat("Printing idxVals[dfrIdxs$pZ50] indices to pMdZ\n");
          # print(idxVals[dfrIdxs$pZ50])
          # pMdZ = vals[idxVals[dfrIdxs$pZ50]];
          # cat("Printing pMdZ\n");
          # print(pMdZ);
          tmPrM2M[array(c(row_idxs,col_idxs),dim=c(nRWs,2))] = prM2M_AscLogistic(vals[idxVals[dfrIdxs$pZ50]],
                                                                                 vals[idxVals[dfrIdxs$pWdZ]],
                                                                                 dfrIdxs$z);
        }
        dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)=="ascnormal");
        if ((nRWs=nrow(dfrIdxs)) > 0){
          #--inputs are pMdZ, pWdZ, z
          col_idxs = dfrIdxs$col_idx;#--from
          row_idxs = dfrIdxs$row_idx;#--to
          tmPrM2M[array(c(row_idxs,col_idxs),dim=c(nRWs,2))] = prM2M_AscNormal(vals[idxVals[dfrIdxs$pMdZ]],
                                                                               vals[idxVals[dfrIdxs$pWdZ]],
                                                                               dfrIdxs$z);
        }
        tmPrM2M[diags] = AD(0);                 #--remove defined self-transition probabilities
        tmPrM2M[diags] = AD(1)-colSums(tmPrM2M);#--assign self-transition probabilities
        lstS[[names(s_)]] = tmPrM2M;
      }#--is_ loop
      lstY[[names(y_)]] = lstS;
    }#--iy_ loop
  } else {
    stop("unrecognized type option for calcGrowth_PrM2M:",info$option);
  }
  if (verbose) cat("Finished calcGrowth_PrM2M().\n")
  return(lstY);
}#--end of function
