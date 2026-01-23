#'
#' @title Calculate mean growth (postmolt size given premolt size)
#' @description Function to calculate mean growth.
#' @param pGrA - mean post-molt size at zGrA
#' @param zGrA - pre-molt size at yielding pGrA as mean post-molt size
#' @param pGrB - mean post-molt size at zGrB
#' @param zGrB - pre-molt size at yielding pGrB as mean post-molt size
#' @param zBs - pre-molt sizes at which to calculate values
#' @return mean growth vector
#'
#' @details The formula used for mean postmolt size (`mnZs`) is
#'
#' mnZs = pGrA\*exp(log(pGrB/pGrA)/log(zGrB/zGrA)\*log(zBs/zGrA));
#'
#' @examplesIf FALSE
#' # example code
#' grMn = grwMeanPostMoltZ1(33,25,150,125,1,seq(25,180,5),seq(25,180,5));
#'
#' @md
#' @export
#'
grwMeanPostMoltZ1<-function(pGrA,zGrA,pGrB,zGrB,zBs){
  mnZs = pGrA*exp(log(pGrB/pGrA)/log(zGrB/zGrA)*log(zBs/zGrA));
  return(mnZs);
}

#'
#' @title Calculate probability of growth (postmolt size given premolt size)
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
#' @details Mean postmolt size uses [grwMeanPostMoltZ1()]:
#'
#' mnZs = grwMeanPostMoltZ1(pGrA,zGrA,pGrB,zGrB,pGrBeta,zBs_from);
#'
#' @examples
#' # example code
#' grM = grwPwrLaw1(33,25,150,125,1,seq(25,180,5),seq(25,180,5));
#'
#' @md
#' @export
#'
grwPwrLaw1<-function(pGrA,zGrA,pGrB,zGrB,pGrBeta,zBs_from,zBs_to,dZ){
  #--mean post-molt sizes
  mnZs = grwMeanPostMoltZ1(pGrA,zGrA,pGrB,zGrB,zBs_from)
  #--mean molt increments
  mnMIs = mnZs - zBs_from;
  #--actual molt increments
  mis = zBs_to - zBs_from;
  #--scaled to gamma distribution alpha parameters
  prGrsU = pgamma(mis+0.5*dZ,shape=mnMIs/pGrBeta,scale=pGrBeta);
  prGrsL = pgamma(mis-0.5*dZ,shape=mnMIs/pGrBeta,scale=pGrBeta);
  prGrs = prGrsU - prGrsL;
  # tbl=tibble::tibble(zBs_from=zBs_from,zBs_to=zBs_to,
  #                    mnZs=mnZs,mnMIs=mnMIs,
  #                    prGrsU=prGrsU,prGrsL=prGrsL,
  #                    prGrs=prGrs);
  # print(tbl,n=Inf)
  return(prGrs);
}

#'
#' @title Calculate mean growth (molt increment or postmolt size given premolt size)
#' @description Function to calculate mean growth.
#' @param pGrA - ln-scale intercept
#' @param pGrB - ln-scale slope
#' @param zBs - pre-molt sizes at which to calculate values
#' @return mean growth vector
#'
#' @details The formula used for mean molt increment/postmolt size (`mnZs`) is
#'
#' mnZs = exp(pGrA+pGrB*log(zBs));
#'
#' @examplesIf FALSE
#' # example code
#'
#' # example code
#' grM = grwMeanPostMoltZ2(33,150,seq(25,180,5));
#'
#' @md
#' @export
#'
grwMeanPostMoltZ2<-function(pGrA,pGrB,zBs){
  mnZs = exp(pGrA+pGrB*log(zBs));
  return(mnZs);
}

#'
#' @title Calculate probability of growth (postmolt size given premolt size)
#' @description Function to calculate probability of postmomlt size given premolt size.
#' @param pGrA - ln-scale intercept
#' @param pGrB - ln-scale slope
#' @param pGrBeta - gamma distribution scale parameter for post-molt variability
#' @param zBs_from - pre-molt sizes at which to calculate values
#' @param zBs_to - post-molt sizes at which to calculate values
#' @return growth matrix
#'
#' @details Mean postmolt size uses [grwMeanPostMoltZ2()]:
#'
#' mnZs = grwMeanPostMoltZ2(pGrA,pGrB,pGrBeta,zBs_from);
#'
#' @examples
#' # example code
#' mtxPrGrw = grwPwrLaw1(33,25,150,125,1,seq(25,180,5),seq(25,180,5));
#'
#' @md
#' @export
#'
grwPwrLaw2<-function(pGrA,pGrB,pGrBeta,zBs_from,zBs_to,dZ){
  #--mean post-molt size
  mnZs = grwMeanPostMoltZ2(pGrA,pGrB,zBs_from);
  #--mean molt increments
  mnMIs = mnZs - zBs_from;
  #--actual molt increments
  mis = zBs_to - zBs_from;
  #--gamma distribution with shape and scale parameters
  prGrs = pgamma(mis+0.5*dZ,shape=mnMIs/pGrBeta,scale=pGrBeta) -
          pgamma(mis-0.5*dZ,shape=mnMIs/pGrBeta,scale=pGrBeta)
  return(prGrs);
}

#' @title Convert parameters/constants values from [grwPwrLaw1()] to [grwPwrLaw2()]
#' @description Function to convert parameters/constants values from [grwPwrLaw1()] to [grwPwrLaw2()]
#' @param pGrA - grwPwrLaw1 parameter pGrA
#' @param pGrB - grwPwrLaw1 parameter pGrB
#' @param zGrA - grwPwrLaw1 constant zGrA
#' @param zGrB - grwPwrLaw1 constant zGrB
#' @param v - vector of grwPwrLaw1 parameters in above order
#' @return vector of grwPwrLaw2 parameters pGrA and pGrB, plus zGrA and zGrB
#' @details Converts parameters/constants for grwPwrLaw1 to those for grwPwrLaw2.
#' @examples
#' par0 = c(32.2, 166.0,  25.0, 125.0); #--pGrA, pGrB, zGrA, zGrB
#' par1 = convertParams_PwrLaw1to2(v=par0);
#' par2 = convertParams_PwrLaw2to1(v=par1);
#' sum(abs(par2-par0));#--should be 0
#'
#' @export
#'
convertParams_PwrLaw1to2<-function(pGrA=NULL,pGrB=NULL,zGrA=NULL,zGrB=NULL,v=NULL){
  if (!is.null(v)){
    pGrA = v[1];
    pGrB = v[2];
    zGrA = v[3];
    zGrB = v[4];
  }
  # f(z) = pGrA*exp(log(pGrB/pGrA)/log(zGrB/zGrA)*log(z/zGrA))
  #      = exp(log(pGrAexp)-log(pGrB/pGrA)/log(zGrB/zGrA)*log(zGrA) +
  #                        +log(pGrB/pGrA)/log(zGrB/zGrA)*log(z))
  pGrAp = log(pGrA)-log(pGrB/pGrA)/log(zGrB/zGrA)*log(zGrA);
  pGrBp = log(pGrB/pGrA)/log(zGrB/zGrA);
  return(c(pGrAp,pGrBp,zGrA,zGrB));
}

#' @title Convert parameters/constants values from [grwPwrLaw2()] to [grwPwrLaw1()]
#' @description Function to convert parameters/constants values from [grwPwrLaw2()] to [grwPwrLaw1()]
#' @param pGrA - grwPwrLaw2 parameter pGrA
#' @param pGrB - grwPwrLaw2 parameter pGrB
#' @param zGrA - grwPwrLaw1 constant zGrA
#' @param zGrB - grwPwrLaw1 constant zGrB
#' @param v - vector of grwPwrLaw2 parameters + grwPwrLaw1 constants in above order
#' @return vector of grwPwrLaw1 parameters pGrA and pGrB, and constants zGrA and zGrB
#' @details Converts parameters/constants for grwPwrLaw2 to those for grwPwrLaw1.
#' @examples
#' par0 = c(0.1919238, 1.0190025, 25.0, 125.0); #--pGrA, pGrB, zGrA, zGrB
#' par1 = convertParams_PwrLaw2to1(v=par0);
#' par2 = convertParams_PwrLaw1to2(v=par1);
#' sum(abs(par2-par0));#--should be 0
#'
#' @export
#'
convertParams_PwrLaw2to1<-function(pGrA=NULL,pGrB=NULL,zGrA=NULL,zGrB=NULL,v=NULL){
  if (!is.null(v)){
    pGrA = v[1];
    pGrB = v[2];
    zGrA = v[3];
    zGrB = v[4];
  }
  # f(z) = exp(pGrA + pGrB*log(z))
  # f(zGrA) = exp(pGrA + pGrB*log(zGrA)) = pGrAp
  # f(zGrB) = exp(pGrA + pGrB*log(zGrB)) = pGrBp
  pGrAp = exp(pGrA+pGrB*log(zGrA));
  pGrBp = exp(pGrA+pGrB*log(zGrB));
  return(c(pGrAp,pGrBp,zGrA,zGrB));
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
#' @return a list of a list of matrices, indexed by year, then season
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
  dZ = unname(dims$zc[2]-dims$zc[1]);#--size bin size
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
      if ("opr_idx" %in% names(dfrUCr))
        if (!is.na(dfrUCr$opr_idx[1])){
          if (dfrUCr$op_type=="additive") {
            p = p + params$pPrGr_OPs[dfrUCr$opr_idx[1]];
          } else {p = p * params$pPrGr_OPs[dfrUCr$opr_idx[1]];}
        }
      if ("dpr_idx" %in% names(dfrUCr))
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
        ####--model categories for y_, s_ (ic_ maps to population state vector index)
        dfrDims  = (dims$dmsYSC |> dplyr::filter(y==y_,s==s_)) |> dplyr::mutate(ic_=dplyr::row_number());
        # idx_base = dfrDims$sparse_idx[1]-1;#--don't need this if using ic_'s as in next
        #
        ####--Determine mapping from pre-molt state to post-molt state
        #####--for each state, z is pre-molt size (also column index).
        #####--then expand to all post-molt sizes (also row indices)
        #####--with same y,s,r,x,m,p (growth would not change any of these categories)
        dfrDimsTF = (dfrDims |> dplyr::rename(z_from=z,col_idx=ic_)) |>
                      dplyr::left_join((dfrDims |> dplyr::rename(z_to=z,row_idx=ic_)),
                                       by=dplyr::join_by(y, s, r, x, m, p),
                                       relationship="many-to-many") |>
                     dplyr::select(y,s,r,x,m,p,z_from,z_to,col_idx,row_idx);

        ####--match functions, parameters info to states that will undergo a growth transition
        #####--NOTE: no negative molt increments allowed, so keep only rows where z_from <= z_to
        dfrIdxsA = dfrDimsTF |>
                     dplyr::inner_join(info$dfrHCs,by = dplyr::join_by(y, s, r, x, m, p, z_from, z_to)) |>
                     dplyr::mutate(z_from=as.numeric(z_from),
                                   z_to  =as.numeric(z_to)) |>
                     dplyr::filter(z_from <= z_to);

        ####--determine (column) indices for states that underwent a growth transition
        #####--need the indices to:
        #####--1. assign a self-transition probability of 1 to states that DID NOT undergo a growth transition
        #####--2. correctly assign the probability of growth into the largest post-molt size bin (an accumulator bin)
        ######--NOTE: this should be calculated in extractParamInfo (it doesn't change for a given y_, s_)
        dfrSUGTs = dfrIdxsA |>
                     dplyr::group_by(y,s,r,x,m,p,z_from,col_idx) |>
                     dplyr::summarize(z_to   =max(z_to),
                                      row_idx=max(row_idx));


        #--function = grwPwrLaw1
        dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)==tolower("grwPwrLaw1"));
        if ((nRWs=nrow(dfrIdxs)) > 0){
          #--parameters are: pGrA,zGrA,pGrB,zGrB,pGrBeta
          col_idxs = dfrIdxs$col_idx;#--from
          row_idxs = dfrIdxs$row_idx;#--to
          tmPrGr[array(c(row_idxs,col_idxs),dim=c(nRWs,2))] = grwPwrLaw1(vals[idxVals[dfrIdxs$pGrA]],
                                                                         vals[idxVals[dfrIdxs$zGrA]],
                                                                         vals[idxVals[dfrIdxs$pGrB]],
                                                                         vals[idxVals[dfrIdxs$zGrB]],
                                                                         vals[idxVals[dfrIdxs$pGrBeta]],
                                                                         dfrIdxs$z_from,dfrIdxs$z_to,dZ);
        }
        #--function = grwPwrLaw2
        dfrIdxs  = dfrIdxsA |> dplyr::filter(tolower(fcn)=="grwPwrLaw2");
        if ((nRWs=nrow(dfrIdxs)) > 0){
          #--parameters are: pGrA,pGrB,pGrBeta
          col_idxs = dfrIdxs$col_idx;#--from
          row_idxs = dfrIdxs$row_idx;#--to
          tmPrGr[array(c(row_idxs,col_idxs),dim=c(nRWs,2))] = grwPwrLaw2(vals[idxVals[dfrIdxs$pGrA]],
                                                                         vals[idxVals[dfrIdxs$pGrB]],
                                                                         vals[idxVals[dfrIdxs$pGrBeta]],
                                                                         dfrIdxs$z_from,dfrIdxs$z_to,dZ);
        }

        #####--for molting categories, calculate transition prob to max post-molt bin
        matAccs = array(c(dfrSUGTs$row_idx,dfrSUGTs$col_idx),dim=c(nrow(dfrSUGTs),2));#--indices of accumulator bins (can be defined in extractParamInfo)
        tmPrGr[matAccs] = AD(0.0);
        csums = colSums(tmPrGr)[dfrSUGTs$col_idx];
        tmPrGr[matAccs] = AD(1.0)-csums;

        ####--for non-molting categories, assign self-transition probabilities of 1
        st_idxs = which(!((1:dims$nCs) %in% dfrSUGTs$col_idx));
        tmPrGr[array(c(st_idxs,st_idxs),dim=c(length(st_idxs),2))] = AD(1.0);
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
