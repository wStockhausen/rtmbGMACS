#--estimate allometry
require(ggplot2);
require(Matrix);
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"testing/newTest_FilesCTL-Allometry");

#--set up model dimensions----
require(rtmbGMACS);
#--source(file.path(dirThs,"r01_SetupDims.R")); #--done in r02_SetupParamVals.R

#--set up model parameters----
##--The following creates a TEMPLATE for specifying initial parameter values and other
##--information. It must be modified from the default values as necessary and copied into
##--the actual CTL file to use. Only modify/run this if creating new models for parameters.
#--source(file.path(dirThs,"r02_SetupParamVals.R"));

#--read in CTL file----
#--source(file.path(dirPrj,"R/readFile_CTL.R"))
##--this sets up dims, parameter values, etc. from the txt file
lstCTL = readFile_CTL(fn=file.path(dirThs,"CTL_Allometry_FormulaApproach.txt"));

#--determine RTMB parameters and initial values for Allometry----
##--get RTMB parameter vectors with initial values----
lstRTMB = rtmbGMACS:::createRTMB_Params(lstCTL$lstAllometryInfo$lstParsInfo);
View(lstRTMB$dfrMM)
fcn_ids = names(lstCTL$lstAllometry$lstFcnInfo$lstFIs);
dfrFI = lstCTL$lstAllometry$lstFcnInfo$lstFIs[[fcn_ids[1]]];
dfrFPI = dfrFI |> dplyr::full_join(lstRTMB$dfrMM)
View(dfrFPI);
lstRTMB$lstRTMB_Xs$fe_pLnA_1 %*% lstRTMB$lstRTMB_Params$fe_pLnA_1;
lstRTMB$lstRTMB_Xs$fe_pB_1 %*% lstRTMB$lstRTMB_Params$fe_pB_1;
lstRTMB$lstRTMB_Xs$fe_pZ0_1 %*% lstRTMB$lstRTMB_Params$fe_pZ0_1;

##--define statistical family to use for fits to allometric data----
family = gaussian(link="log");

##--define dataframe for "standard" weight-at-size prediction
dfrPrd = mfYXM |> dplyr::cross_join(tibble::tibble(z=seq(25,180,5)));

##--collect allometric information----
lstAllom = list(dataDFR=dfrData,
                dfrPrd=dfrPrd,
                lstRTMB=lstRTMB
                );

#--create input list----
###--list will be expanded in complete model objective function with info for other processes
inputs = list(lstAllom=lstAllom
             );

params = lstRTMB$lstRTMB_Params;
re_nms = lstCTL$reParamNames;
map    = lstRTMB$lstRTMB_Map;

rm(lstAllom,lstRTMB);

#--source objective function----
source(file.path(dirThs,"objfun_EstimateAllometry.R"));


#--run objective function in R environment for testing----
if (FALSE) {
  testing=TRUE;

  ##--run the function
  objfun<-objfun_EstimateAllometry(params);

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
  objfun<-MakeADFun(objfun_EstimateAllometry,
                    params,
                    random=re_nms,
                    map=map,
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

