#--estimate allometry
require(ggplot2);
dirThs = dirname(rstudioapi::getSourceEditorContext()$path);

#--get allometric data----
dfrZW  = wtsUtilities::getObj(file.path(dirThs,"rda_AllometryData.RData")) |>
           dplyr::rename(obs=w) |> tibble::rownames_to_column(var="obs_id") |>
           dplyr::mutate(dplyr::across(dplyr::everything(),as.character));
plt = ggplot(dfrZW,aes(x=as.numeric(z),y=as.numeric(obs),colour=x)) + geom_point();
print(plt);
##--remove some clearly bad values
dfrZWp = dfrZW |>
          dplyr::filter(!((as.numeric(z)>125)&(as.numeric(obs)<100))) |>
          dplyr::filter(!((as.numeric(z)<100)&(as.numeric(obs)>400))) |>
          dplyr::filter(!((as.numeric(z)< 40)&(as.numeric(obs)> 75))) |>
          dplyr::filter(as.numeric(obs)>0);
plt = ggplot(dfrZWp,aes(x=as.numeric(z),y=as.numeric(obs),colour=x)) + geom_point() +
        coord_cartesian(xlim=c(0,60),ylim=c(0,100))
print(plt);
dfrZW = dfrZWp;

#--set up model dimensions----
require(rtmbGMACS);
source(file.path(dirThs,"r01_SetupDims.R"));

#--read in CTL file----
str = readr::read_lines(file.path(dirThs,"CTL_Allometry1.txt")) |>
      paste("\n",collapse=""); #--need this as a single string

#--process CTL file----
strv = str |> splitText() |> extractLines("PROCESS_ALLOMETRY","END_ALLOMETRY");
strv |> extractLines("CODE","END_CODE") |> evalTextAsCode(frame=.GlobalEnv);#--change to 0 in objfun

##--process function information----
dfrFcns = strv |> extractLines("FUNCTIONS","END_FUNCTIONS")   |> evalTextAsDataframe();
###--expand function info to relevant model dimensions----
lstFcnsExp = list();
for (row in 1:nrow(dfrFcns)){
  rw = dfrFcns[row,];
  lstFcnsExp[[row]] = get(rw$frame) |> dplyr::mutate(fcn_idx=rw$fcn_idx,`function`=rw$`function`,params=list(stringr::str_split_1(rw$params,",")));
}
dfrFcnsExp = dplyr::bind_rows(lstFcnsExp);
rm(row,rw,lstFcnsExp);

##--process parameters equation information
dfrPEQs = strv |> extractLines("PARAM_EQS","END_PARAM_EQS")   |> evalTextAsDataframe();
par_ids = dfrPEQs$par_id;
###--fixed effects parameter values
strFEVs = strv |> extractLines("FE_VALS","END_FE_VALS");
lstFEVs = list();
if (length(strFEVs)>0){
  for (par_id in par_ids){
    strp = strFEVs |> extractLines(par_id,paste0("END_",par_id));
    if (length(strp)>0)
      lstFEVs[[par_id]]  = strp |> evalTextAsDataframe() |> dplyr::mutate(par_id=par_id,.before=1) |> expand_mirrors();
    rm(strp);
  }
  rm(par_id);
}#--FEVs
###--environmental covariates parameter values
strECVs = strv |> extractLines("EC_VALS","END_EC_VALS");
lstECVs = list();
if (length(strECVs)>0){
  for (par_id in par_ids){
    strp = strECVs |> extractLines(par_id,paste0("END_",par_id));
    if (length(strp)>0)
      lstECVs[[par_id]]  = strp |> evalTextAsDataframe() |> dplyr::mutate(par_id=par_id,.before=1) |> expand_mirrors();
    rm(strp);
  }
  rm(par_id);
}#--ECVs
###--random effects parameter values
strREVs = strv |> extractLines("RE_VALS","END_RE_VALS");
lstREVs = list();
if (length(strREVs)>0){
  for (par_id in par_ids){
    strp = strREVs |> extractLines(par_id,paste0("END_",par_id));
    if (length(strp)>0)
      lstREVs[[par_id]]  = strp |> evalTextAsDataframe() |> dplyr::mutate(par_id=par_id,.before=1) |> expand_mirrors();
    rm(strp);
  }
  rm(par_id);
}#--REVs

##--determine model matrices for conversion from actual RTMB parameters to model parameters
tbl = dfrFcns |> dplyr::full_join(dfrPEQs,by="fcn_idx");
idxFEs = idxECs = idxREs = 0;
lstModMtx = list();
for (par_idx_ in tbl$par_idx){
  #--testing: par_idx_ = 3;
  rw = tbl |> dplyr::filter(par_idx==par_idx_);
  #--fixed effects
  lstModMtxFEs = calcModelMatrixFEs(txt=rw$feEQs,dfrMF=get(rw$frame),ctrs=get(rw$feContrasts));
  idxFEs = idxFEs + lstModMtxFEs$npars;
  lstModMtxFEs$idx_end = idxFEs;
  #--env. covariates
  lstModMtxECs = calcModelMatrixFEs(txt=rw$ecEQs,dfrMF=get(rw$frame),ctrs=NULL,verbose=TRUE);
  idxECs = idxECs + lstModMtxECs$npars;
  lstModMtxECs$idx_end = idxECs;
  #--random effects
  lstModMtxREs = calcModelMatrixREs(rw$reEQs,dfrMF=get(rw$frame),cov_type=rw$reCovType);
  idxREs = idxREs + lstModMtxREs$npars;
  lstModMtxREs$idx_end = idxREs;
  #--combine lists into one object
  nm=paste0("allom_",rw$par_id,"-",rw$par_idx);
  lstModMtx[[nm]] = list(par_nm=nm,par_idx=rw$par_idx,par_id=rw$par_id,links=getLinkFcn(rw$link_fcn),
                         FEs=lstModMtxFEs,ECs=lstModMtxECs,REs=lstModMtxREs);
  rm(rw,lstModMtxFEs,lstModMtxECs,lstModMtxREs);
}
rm(par_idx_,nm);

##--map functions/parameters to data----
dfrZW1 = dfrZW |>
           dplyr::filter(dplyr::between(as.numeric(y),2015,2024)) |>  #--keep only data w/in model time frame5
           dplyr::mutate(m=ifelse(m=="all",NA,m),     #--convert maturity state "all" to NA
                         p=ifelse(x=="female",NA,p));#--convert post-molt age for females to NA
dfrZWp = dfrZW1  |>
           dplyr::left_join(dfrFcnsExp) |>
           dplyr::arrange(obs_id);
lstZWp = list();
for (nm in names(lstModMtx)){
  #--testing: nm = names(lstModMtx)[1];
  lstModMtxPar = lstModMtx[[nm]];
  if (lstModMtxPar$FEs$npars>0) {
    dfrZWppFEs = dfrZWp |> dplyr::mutate(par_nm=nm,par_id=lstModMtxPar$par_id,par_idx=lstModMtxPar$par_idx,type="FEs") |>
                   dplyr::inner_join(lstModMtxPar$FEs$dfrMF);
  } else {dfrZWppFEs = NULL;}
  if (lstModMtxPar$ECs$npars>0) {
    dfrZWppECs = dfrZWp |> dplyr::mutate(par_nm=nm,par_id=lstModMtxPar$par_id,par_idx=lstModMtxPar$par_idx,type="ECs") |>
                   dplyr::inner_join(lstModMtxPar$ECs$dfrMF);
  } else {dfrZWppECs = NULL;}
  if (lstModMtxPar$REs$npars>0) {
    dfrZWppREs = dfrZWp |> dplyr::mutate(par_nm=nm,par_id=lstModMtxPar$par_id,par_idx=lstModMtxPar$par_idx,type="REs") |>
                   dplyr::inner_join(lstModMtxPar$REs$dfrMF);
  } else {dfrZWppREs = NULL;}
  lstZWp[[nm]] = dplyr::bind_rows(dfrZWppFEs,dfrZWppECs,dfrZWppREs);# |> tidyr::pivot_wider(names_from="type",values_from="row_idx");
  rm(lstModMtxPar,dfrZWppFEs,dfrZWppECs,dfrZWppREs);
}
dfrZWpp = dplyr::bind_rows(lstZWp) |>
            tidyr::pivot_wider(names_from="type",values_from="row_idx") |>
            dplyr::arrange(obs_id);
if (!("FEs" %in% names(dfrZWpp))) dfrZWpp$FEs = NA;
if (!("ECs" %in% names(dfrZWpp))) dfrZWpp$ECs = NA;
if (!("REs" %in% names(dfrZWpp))) dfrZWpp$REs = NA;
#dfrZWpp  = dfrZWpp |> dplyr::rowwise() |> dplyr::mutate(row_idxs = list(c(FEs,ECs,REs)));#--NOTE: might want to drop row_idxs?
udfrZWpp = dfrZWpp |> dplyr::distinct(fcn_idx,par_idx,FEs,ECs,REs);#--unique function/model parameter combinations
rm(nm);

##--define statistical family to use for fits to allometric data----
family = gaussian(link="log");

##--collect allometric information----
lstAllom = list(dfrData=dfrZW, mapData=dfrZWpp, uMapData=udfrZWpp, family=family, lstModMtx=lstModMtx, nparsFEs=idxFEs, nparsECs=idxECs, nparsREs=idxREs);

#--determine RTMB parameters and initial values----
#--create map for parameter phasing, mirroring
##--get priors info for model parameters
##--determine initial values for "actual" RTMB parameters
###--using svd inversion and link function to determine initial values
###--of "actual" RTMB parameters from input values for model parameters
lst_pinfo  = list();
lst_map    = list();
lst_priors = list();
params     = list();
##--process fixed effects info----
if (length(lstFEVs)>0){
  for (par_id_ in names(lstFEVs)){
    ##--testing: par_id_ = names(lstFEVs)[1];
    dfrIVs   = lstFEVs[[par_id_]];               #--dataframe with model parameter info
    for (par_idx_ in unique(dfrIVs$par_idx)){
      #--testing: par_idx_ = unique(dfrIVs$par_idx)[1];
      nm = paste0("allom_",par_id_,"-",par_idx_);
      dfrIVsp = dfrIVs |> dplyr::filter(par_idx_==par_idx) |> dplyr::mutate(par_nm=paste0(nm,"_FEs"));
      ##--create "actual" RTMB parameter vector corresponding to model parameter `par_id_`
      nmp = paste0(nm,"_FEs");
      lst_pinfo[[nmp]]  = getParametersInitialValues(lst=lstAllom$lstModMtx[[nm]]$FEs,dfrIVsp,links=lstAllom$lstModMtx[[nm]]$links);
      params[[nmp]]     = lst_pinfo[[nmp]]$vP;
      lst_map[[nmp]]    = createParamsMap(dfrIVsp);       #--create RTMB "map" for model parameter vector
      lst_priors[[nmp]] = getParametersPriorInfo(dfrIVsp);#--get prior info for parameters
    } #--loop: par_idx_ in unique(dfrIVs$par_idx)
  } #--loop: par_id_ in names(lstFEVs)
} #--if: length(lstFEVs)>0

##--process environmental covariates info----
if (length(lstECVs)>0){
  for (par_id_ in names(lstECVs)){ #--loop over model parameter names
    ##--testing: par_id_ = names(lstECVs)[1];
    dfrIVs   = lstECVs[[par_id_]];  #--dataframe with model parameter info
    for (par_idx_ in unique(dfrIVs$par_idx)){ #--loop over "indices" associated with the model parameter
      #--testing: par_idx_ = unique(dfrIVs$par_idx)[1];
      nm = paste0("allom_",par_id_,"-",par_idx_);#--RTMB parameter base name
      dfrIVsp = dfrIVs |> dplyr::filter(par_idx_==par_idx) |> dplyr::mutate(par_nm=paste0(nm,"_ECs"));
      ##--create "actual" RTMB parameter vector corresponding to model parameter `par_id_`
      nmp = paste0(nm,"_ECs");
      lst_pinfo[[nmp]]  = getParametersInitialValues(lst=lstAllom$lstModMtx[[nm]]$ECs,dfrIVsp=dfrIVsp,links=lstAllom$lstModMtx[[nm]]$links);
      params[[nmp]]     = (lst_pinfo[[nmp]])$vP;
      lst_map[[nmp]]    = createParamsMap(dfrIVsp);       #--create RTMB "map" for model parameter vector
      lst_priors[[nmp]] = getParametersPriorInfo(dfrIVsp);#--get prior info for parameters
    }#--par_idx_ in unique(dfrIVs$par_idx) loop
  }#--par_id_ in names(lstECVs)
}#--end processing lstECVs

##--process random effects info  TODO!!----

##--finish constructing lstAllom----
lstAllom[["lst_pinfo"]]  = lst_pinfo;
lstAllom[["lst_map"]]    = lst_map;
lstAllom[["lst_priors"]] = lst_priors;

dfrZWppp = dfrZWpp |>
             dplyr::rowwise() |>
             dplyr::mutate(par_info=list(list(par_id=par_id,par_idx=par_idx,FEs=FEs,ECs=ECs,REs=REs))) |>
             dplyr::ungroup() |>
             dplyr::select(!c(par_nm,par_idx,FEs,ECs,REs)) |>
             tidyr::pivot_wider(names_from="par_id",values_from=par_info) |>
             dplyr::mutate(pars=list(list()),
                           prd=NA_real_,
                           nll=NA_real_);

inputs = list(lstAllom=lstAllom);

#--source objective function----
source(file.path(dirThs,"objfun_EstimateAllometry.R"));

plt = ggplot(dfrZWpp |> dplyr::select(z,obs,x),aes(x=as.numeric(z),y=as.numeric(obs),colour=x)) + geom_point();
print(plt);
plt + scale_y_log10() + scale_x_log10();

#--run objective function in R environment for testing----
if (FALSE) {
  testing=TRUE;
  objfun<-objfun_EstimateAllometry(params);
}

#--run in RTMB environment----
if (FALSE){
  require(RTMB);
  testing=FALSE;
  objfun<-MakeADFun(objfun_EstimateAllometry,
                    params,
                    random=NULL,
                    map=list(),
                    ADreport=FALSE,
                    ridge.correct=FALSE, #--Note: apparently EXPERIMENTAL
                    silent=FALSE);
  opt <- nlminb(objfun$par, objfun$fn, objfun$gr);
  rep <- sdreport(objfun);
  summary(rep);
  summary(rep,"fixed");
  summary(rep,"random");
  summary(rep,"report");

  rep = objfun$report(objfun$env$last.par.best);
  dfrZWs = rep$dfrZWs;
  plt = ggplot() +
          geom_point(data=dfrZWs,mapping=aes(x=obs,y=prd,colour=x,fill=x)) +
          geom_abline(slope=1);
  print(plt);
  plt = ggplot() +
          geom_point(data=dfrZWs,mapping=aes(x=z,y=obs,colour=x,fill=x));
  print(plt);
}

