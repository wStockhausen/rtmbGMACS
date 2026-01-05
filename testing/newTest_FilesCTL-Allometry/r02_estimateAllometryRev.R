#--estimate allometry
require(ggplot2);
dirPrj = rstudioapi::getActiveProject();
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
rm(dfrZWp);

#--set up model dimensions----
require(rtmbGMACS);
source(file.path(dirThs,"r01_SetupDims.R"));

#--read in CTL file----
source(file.path(dirPrj,"R/readFile_CTLrev.R"))
lstCTL = readFile_CTLrev(fn=file.path(dirThs,"CTL_Allometry3rev.txt"));

##--determine model matrices for conversion from actual RTMB parameters to model parameters
source(file.path(dirPrj,"R/determineModelMatricesRev.R"))
lstModMtx = determineModelMatricesRev(lstCTL$dfrPEQs);

##--determine map from model parameter indices to functions
source(file.path(dirPrj,"R/mapParameterIndicesToFunctions.R"))
dfrMapPIsToFcns = mapParameterIndicesToFunctions(lstCTL,lstModMtx);


#--determine RTMB parameters and initial values----
#--create map for parameter phasing, mirroring
##--get priors info for model parameters
##--determine initial values for "actual" RTMB parameters
source(file.path(dirPrj,"R/determineInitialValuesRev.R"))
source(file.path(dirPrj,"R/MiscFunctions_Parameters.R")) #--for `createParamsMapRev`
lstInfoIVs = determineInitialValuesRev(lstCTL,lstModMtx);

##--map functions/parameters to data----
dfrZW1 = dfrZW |>
           dplyr::filter(dplyr::between(as.numeric(y),2015,2024)) |>  #--keep only data w/in model time frames
           dplyr::mutate(m=ifelse(m=="all",NA,m),     #--convert maturity state "all" to NA
                         p=ifelse(x=="female",NA,p)); #--convert post-molt age for females to NA
dfrZWppp = dfrZW1  |>
             dplyr::left_join(dfrMapPIsToFcns) |>
             dplyr::arrange(obs_id);

##--define statistical family to use for fits to allometric data----
family = gaussian(link="log");

##--define dataframe for "standard" weight-at-size prediction
dfrPrd = dfrMapPIsToFcns |> dplyr::cross_join(tibble::tibble(z=seq(25,180,5)));

##--collect allometric information----
lstAllom = list(dataDFR=dfrZWppp,
                dfrPrd=dfrPrd,
                dfrMapPIsToFcns=dfrMapPIsToFcns,
                family=family);
rm(dfrZWppp,dfrMapPIsToFcns,family);

#--create input list----
###--list will be expanded in complete model objective function with info for other processes
inputs = list(lstAllom=lstAllom,
                lstModMtx=lstModMtx,
                lst_pinfo=lstInfoIVs$lst_pinfo,
                lst_map=lstInfoIVs$lst_map,
                lst_priors=lstInfoIVs$lst_priors);

params = lstInfoIVs$params;

rm(lstAllom,lstModMtx,lstInfoIVs);

#--source objective function----
source(file.path(dirThs,"objfun_EstimateAllometryRev.R"));

plt = ggplot(inputs$lstAllom$dataDFR |> dplyr::select(z,obs,x) |> dplyr::mutate(z=as.numeric(z),obs=as.numeric(obs)),
             aes(x=z,y=obs,colour=x)) + geom_point();
print(plt);
plt + scale_y_log10() + scale_x_log10();

#--run objective function in R environment for testing----
if (FALSE) {
  testing=TRUE;
  objfun<-objfun_EstimateAllometryRev(params);
  plt = ggplot(objfun$dfrZWs |> dplyr::select(z,obs,prd,x),aes(x=z,y=obs,colour=x)) +
          geom_point() + geom_point(aes(y=prd),colour="black",size=1) + facet_grid(x~.);
  print(plt);
  plt + scale_y_log10() + scale_x_log10();
  plt = ggplot(objfun$dfrZWs |> dplyr::select(z,obs,prd,x) |> dplyr::mutate(residual=log(obs/prd)),aes(x=z,y=residual,colour=x)) +
          geom_point() + facet_grid(x~.);
  print(plt);
}

#--run in RTMB environment----
if (TRUE){
  require(RTMB);
  testing=FALSE;
  objfun<-MakeADFun(objfun_EstimateAllometryRev,
                    params,
                    random=NULL,
                    map=inputs$lst_map["allom_pZ0_1_FEs"],
                    hessian=TRUE,
                    LaplaceNonZeroGradient=TRUE,
                    ADreport=FALSE,
                    ridge.correct=FALSE, #--Note: apparently EXPERIMENTAL
                    silent=FALSE);
  opt <- nlminb(objfun$par, objfun$fn, objfun$gr,objfun$he);
  sdrep <- sdreport(objfun);
  print(sdrep);
  mPrd = summary(sdrep,"report");
  dfrPrdRes = inputs$lstAllom$dfrPrd |>
                dplyr::bind_cols(tibble::as_tibble(mPrd[rownames(mPrd)=="prd_allom",])) |>
                dplyr::mutate(ymin=Estimate-1.96*`Std. Error`,ymax=Estimate+1.96*`Std. Error`);
  ggplot(dfrPrdRes |> dplyr::filter(x=="male",yn==2015),aes(x=as.numeric(z),y=Estimate,colour=m,group=as.factor(yn))) +
    geom_line() + geom_ribbon(aes(ymin=ymin,ymax=ymax),alpha=0.3) +
    scale_y_log10() + labs(x="size (mm CW)",y="predicted weight (gm)");
  ggplot(dfrPrdRes |> dplyr::filter(x=="female",yn==2015),aes(x=as.numeric(z),y=Estimate,group=as.factor(yn))) +
    geom_line() + geom_ribbon(aes(ymin=ymin,ymax=ymax),alpha=0.3) +
    facet_grid(~m) +
    scale_y_log10() + labs(x="size (mm CW)",y="predicted weight (gm)");

  par_best = objfun$env$last.par.best;
  objfun$fn(par_best);
  objfun$gr(par_best);
  objfun$he(par_best);
  rep = objfun$report(par_best);
  dfrZWs = rep$dfrZWs;
  plt = ggplot() +
          geom_point(data=dfrZWs,mapping=aes(x=obs,y=prd,colour=x,fill=x)) +
          geom_abline(slope=1);
  print(plt);
  plt = ggplot() +
          geom_point(data=dfrZWs,mapping=aes(x=z,y=nll,colour=x,fill=x)) + facet_grid(x~.);
  print(plt);
  plt = ggplot(dfrZWs |> dplyr::select(z,obs,prd,x) |> dplyr::mutate(residual=log(obs/prd)),aes(x=z,y=residual,colour=x)) +
          geom_point() + facet_grid(x~.);
  print(plt);
}

