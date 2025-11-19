
getParamInfo<-function(strv,process_name="allom"){
  lstPI = list(); #--list for param info output

  #--extract formula-type equations defining parameter structures
  dfrPEQs = strv |> extractLines("PARAM_EQS","END_PARAM_EQS") |> evalTextAsDataframe();
  par_ids = dfrPEQs$par_id;

  #--extract parameter values
  ##--fixed effects parameter values
  strFEVs = strv |> extractLines("FE_VALS","END_FE_VALS");
  lstFEVs = list();
  if (length(strFEVs)>0){
    for (par_id in par_ids){
      strp = strFEVs |> extractLines(par_id,paste0("END_",par_id));
      lstFEVs[[par_id]]  = strp |> evalTextAsDataframe() |> dplyr::mutate(par_id=par_id,.before=1) |> expand_mirrors();
    }
  }#--FEVs
  ##--environmental covariates parameter values
  strECVs = strv |> extractLines("EC_VALS","END_EC_VALS");
  lstECVs = list();
  if (length(strECVs)>0){
    for (par_id in par_ids){
      strp = strECVs |> extractLines(par_id,paste0("END_",par_id));
      lstECVs[[par_id]]  = strp |> evalTextAsDataframe() |> dplyr::mutate(par_id=par_id,.before=1) |> expand_mirrors();
    }
  }#--ECVs
  ##--random effects parameter values
  strREVs = strv |> extractLines("RE_VALS","END_RE_VALS");
  lstREVs = list();
  if (length(strREVs)>0){
    for (par_id in par_ids){
      strp = strREVs |> extractLines(par_id,paste0("END_",par_id));
      lstREVs[[par_id]]  = strp |> evalTextAsDataframe() |> dplyr::mutate(par_id=par_id,.before=1) |> expand_mirrors();
    }
  }#--REVs

  #--determine model matrices
  tbl = dfrFcns |> dplyr::full_join(dfrPEQs,by="fcn_idx");
  idxFEs = idxECs = idxREs = 0;
  lstModMtx = list();
  for (par_idx_ in tbl$par_idx){
    #--testing: par_idx_ = 1;
    rw = tbl |> dplyr::filter(par_idx==par_idx_);
    #--fixed effects
    lstModMtxFEs = calcModelMatrixFEs(rw$feEQs,dfrMF=get(rw$frame),ctrs=get(rw$feCons));
    idxFEs = idxFEs + lstModMtxFEs$npars;
    lstModMtxFEs$idx_end = idxFEs;
    #--env. covariates
    lstModMtxECs = calcModelMatrixFEs(rw$ecEQs,dfrMF=get(rw$frame));
    idxECs = idxECs + lstModMtxECs$npars;
    lstModMtxECs$idx_end = idxECs;
    #--random effects
    lstModMtxREs = calcModelMatrixREs(rw$reEQs,dfrMF=get(rw$frame),cov_type=rw$reCovType);
    idxREs = idxREs + lstModMtxREs$npars;
    lstModMtxREs$idx_end = idxREs;
    #--combine lists into one object
    lstModMtx[[rw$par_id]] = list(FEs=lstModMtxFEs,ECs=lstModMtxECs,REs=lstModMtxREs);
  }
  lstAllom = list(lstModMtx=lstModMtx, nparsFEs=idxFEs, nparsECs=idxECs, nparsREs=idxREs);

  #--create parameter vectors and parameter map from parameter values dataframes
  lst_params = list();
  lst_map    = list();
  vec_REs    = character(0);
  if (length(lstFEVs)>0){
    for (par_id_ in names(lstFEVs)){
      #--testing: par_id_ = names(lstFEVs)[1];
      parname = paste0(process_name,"_FE_",par_id_);
      dfr = lstFEVs[[par_id_]];
      lst_params[[parname]] = dfr$IV;#--initial values
      lst_map    = c(lst_map,createParamsMap(dfr));
    }
  }
  if (length(lstECVs)>0){
    for (par_id_ in names(lstECVs)){
      #--testing: par_id_ = names(lstECVs)[1];
      parname = paste0(process_name,"_EC_",par_id_);
      dfr = lstECVs[[par_id_]];
      lst_map    = c(lst_map,createParamsMap(dfr));
      lst_map[[parname]]    = createParamsMap(dfr);
    }
  }
  if (length(lstREVs)>0){
    for (par_id_ in names(lstREVs)){
      #--testing: par_id_ = names(lstREVs)[1];
      parname = paste0(process_name,"_RE_",par_id_);
      dfr = lstREVs[[par_id_]];
      lst_params[[parname]] = dfr$IV;#--initial values
      lst_map    = c(lst_map,createParamsMap(dfr));
      vec_REs = c(vec_REs,parname);
    }
  }

  #--create output
  return(lstPI);
}

#' @title Get canonical names for a parameter values dataframe
#' @export
getParamValsDataframe_colnames<-function(){
  all_cols = c("pv_idx",
               "par_idx",
               "mirror","phase","units",
               "IV","LB","UB","prior","p1","p2")
  return(all_cols);
}
#' @title Expand parameter values dataframe for missing columns
#' @export
expandParamDataframe_columns<-function(dfr){
  all_cols = getParamValsDataframe_colnames();
  colnms = all_cols[!(all_cols %in% names(dfr))];
  for (colnm in colnms) dfr = dfr |> tibble::add_column("{colnm}":=NA)
  return(dfr);
}

#' @title Expand dataframe for mirrored parameters
#' @export
expand_mirrors<-function(dfr){
  if (!is.null(dfr$mirror)){
    nm = which(tolower(names(dfr))=="mirror");
    nc  = ncol(dfr);
    tst = dfr$mirror>0;
    for (i in 1:length(tst)){
      if (tst[i]){
        im = dfr$mirror;#--index to mirror
        dfr[i,(nm+1):nc] = dfr[im,(nm+1):nc];
      }
    }
  }
  return(dfr);
}

#' @title Create TMB `map` for parameters from a parameter values dataframe
#' @export
createParamsMap<-function(dfr){
  dfrUPs = dfr |> dplyr::distinct(par_id,par_idx);
  map = list();
  for (r in 1:nrow(dfrUPs)){
    rw = dfrUPs[r,];
    nm = paste0(rw$par_id,"-",rw$par_idx);
    dfrp = dfr |> dplyr::inner_join(rw);
    fac = factor(as.character(dfrp$pv_idx));
    fac = ifelse(dfr$mirror>0,dfrp$mirror,fac);
    fac = ifelse(dfrp$phase<1, NA,fac);
    map[[nm]] = fac;
  }
  return(map);
}
#createParamsMap(dfr);

#' @title Create dataframe with priors info for parameters from a parameter values dataframe
#' @export
getParametersPriorInfo<-function(dfr){
  prior = ifelse(is.null(dfr$prior),rep(NA,nrow(dfr)),"zero");
  p1    = ifelse(is.null(dfr$p1),   rep(NA,nrow(dfr)),as.numeric(dfr$p1));
  p2    = ifelse(is.null(dfr$p2),   rep(NA,nrow(dfr)),as.numeric(dfr$p2));
  dfrp = tibble::tibble(par_id=dfr$par_id,par_idx=dfr$par_idx,prior=prior,p1=p1,p2=p2);
  return(dfrp);
}
#getParametersPriorInfo(dfr);

#'
#' @title Get initial values for a RTMB parameter
#' @description Function to get initial values for a RTMB parameter.
#' @param lst -
#' @param dfrIVsp - dataframe with initial model parameter values
#' @param links - links list (from [getLinkFcn()]) for transformation from/to model to/from RTMB parameter scales, or NULL
#' @param verbose - flag to print diagnostic info
#' @return vector with initial values for RTMB parameter
#' @details TODO!
#' @import glue
#' @importFrom rlang data_sym
#' @export
getParametersInitialValues<-function(lst,dfrIVsp,links=NULL,verbose=FALSE){
  #--get factor names from model matrix/frame
  facnms = lst$factors; #--should this include covariates?
  #--get variable (covariate) names from model matrix/frame
  varnms = lst$variables;
  if (verbose) {
    cat("facnms =",facnms,"\n");
    cat("varnms =",varnms,"\n");
  }
  #--get initial values
  dfrIVspp = dfrIVsp[,c(facnms,"IV")];            #--initial values for model parameter vector
  #--convert NA's to "all"s, then expand initial values as necessary for "all"s
  for (facnm in facnms){
    #--testing: facnm = facnms[1];
    lstLevs = list();
    lstLevs[[facnm]] = (lst$dfrMF |> dplyr::select(dplyr::any_of(facnm)) |> dplyr::distinct())[[facnm]];
    facnm_sym = rlang::data_sym(facnm);
    dfrIVspp = dfrIVspp |> dplyr::mutate({{facnm}}:=ifelse(is.na(!!facnm_sym),"all",!!facnm_sym)) |> alls_ExpandInDataframe(lstLevs);
  }
  vIVs = dfrIVspp$IV;

  ##--calculate reduced model matrix corresponding to model parameter factor levels
  colnms=colnames(lst$mtxRedMM);
  dfrRedT = lst$dfrMF |>
              #dplyr::bind_cols(as.matrix(lst$mtxRedMM)) |>
              dplyr::inner_join(lst$dfrModMtx,by=c(facnms,varnms)) |>
              dplyr::group_by(dplyr::pick(dplyr::any_of(facnms))) |>
              dplyr::summarize(dplyr::across(dplyr::all_of(c(colnms,varnms)),mean)) |> dplyr::ungroup();
  mtxRedT = as.matrix(dfrRedT |> dplyr::select(!dplyr::any_of(c(facnms))));
  ##--get SVD of the reduced model matrix
  lstSVD = svd(mtxRedT);
  ##--get link and inverse link functions defining parameter transformation
  links = lst$links;
  if (is.null(links)) links = getLinkFcn("ident");
  ##--calculate initial values of "actual" RTMB parameters
  vP = as.vector((lstSVD$v %*% diag(1/lstSVD$d,nrow=length(lstSVD$d)) %*% t(lstSVD$u)) %*% links$link(vIVs));
  ##--recalculate initial values of model parameters as check
  vpIVs = as.vector(links$link_inv(mtxRedT %*% vP));
  if (sum(abs(vIVs-vpIVs))>0.0001){
    str = paste0(nm,": ",lst$f,"\n",
                 "  actual initial values:",vP,"\n",
                 "original initial values:",vIVs,"\n",
                 " checked initial values:",vpIVs,"\n\n");
    warning(str);
  }
  return(list(vP=vP,mtxRedT=mtxRedT,links=links));
}

