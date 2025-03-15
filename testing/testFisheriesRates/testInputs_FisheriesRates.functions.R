#--test "function" fishery rates specifications for the gmacs objective function
require(RTMB);
require(ggplot2);
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
setwd(dirThs);
if (FALSE){
  require(rtmbGMACS);
} else {
  source(file.path(dirPrj,"R","DimensionsFunctions.R"))
  source(file.path(dirPrj,"R","DimensionsUtilities.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Alls.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Dataframe.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Text.R"))
  source(file.path(dirPrj,"R","MiscFunctions_Transforms.R"))
  source(file.path(dirPrj,"R","MiscFunctions.R"))
  source(file.path(dirPrj,"R","readParamInfoSectionType1.R"))
  source(file.path(dirPrj,"R","readParamInfo_FisheriesRates.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_FisheriesRates.R"))
  source(file.path(dirPrj,"R","Functions_Catchability.R"))
  source(file.path(dirPrj,"R","calcFisheriesRates.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "function";
if (type=="function"){
  dims = setupModelDims(zcs=seq(25.5,84.5,5));
  params = list();
  map = list();
  ###--selectivity functions with function type input----
  source(file.path(dirPrj,"R","readParamInfo_Selectivity.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Selectivity.R"))
  source(file.path(dirPrj,"R","Functions_Selectivity.R"))
  source(file.path(dirPrj,"R","calcSelectivity.R"))
  conn   = file.path(dirPrj,"testing/testSelectivity/inputSpecs_Selectivity.function.txt");
  res    = readParamInfo_Selectivity(conn,FALSE);
  lstSel = extractParamInfo_Selectivity(res,dims,FALSE);
  params[["pSel_MPs"]]=lstSel$MPs$params;
  map = addList(map,lstSel$map);
  if (!is.null(lstSel$OPs$params)) params[["pSel_OPs"]]=lstSel$OPs$params;
  if (!is.null(lstSel$DPs$params)) params[["pSel_DPs"]]=lstSel$DPs$params;
  if (!is.null(lstSel$REs$params)) params[["pSel_REs"]]=lstSel$REs$params;

  ###--survey catchability inputs with function format----
  conn   = file.path(dirPrj,"testing/testFisheriesRates/inputSpecs_FisheriesRates.function.txt");
  res    = readParamInfo_FisheriesRates(conn,FALSE);
  lstSrv = extractParamInfo_FisheriesRates(res,dims,FALSE);
  params[["pSrv_MPs"]]=lstSrv$MPs$params;
  map = addList(map,lstSrv$map);
  if (!is.null(lstSrv$OPs$params)) params[["pSrv_OPs"]]=lstSrv$OPs$params;
  if (!is.null(lstSrv$DPs$params)) params[["pSrv_DPs"]]=lstSrv$DPs$params;
  if (!is.null(lstSrv$REs$params)) params[["pSrv_REs"]]=lstSrv$REs$params;
}
inputs        = list();
inputs$dims   = dims;
inputs$lstSel = lstSel;#--add lstSel to inputs
inputs$lstSrv = lstSrv;#--add lstSrv to inputs

#--test calcFisheriesRates function----
source(file.path(dirPrj,"R/calcSelectivity.R"));
lstSelVals = calcSelectivity(inputs$dims,inputs$lstSel,params,FALSE);
source(file.path(dirPrj,"R/calcFisheriesRates.R"));
lstSrvVals = calcFisheriesRates(inputs$dims,inputs$lstSrv,params,lstSelVals,TRUE);

#--plot results----
plotSrvVals<-function(if_=1,iy_=1,is_=1){
  dfrSrv = dims$dmsC |>
            dplyr::mutate(value=(lstSrvVals[[if_]])[iy_,is_,],
                          z=as.numeric(as.character(z)));
  ggplot(dfrSrv, aes(x=as.numeric(z),y=value)) +
    geom_line() + geom_hline(yintercept=0,linetype=3,colour="grey") +
    facet_grid(m+p~x) +
    labs(x="size",y="survey catchability") +
    wtsPlots::getStdTheme();
}
plotSrvVals(1,1,1)
plotSrvVals(2,1,1)
plotSrvVals(1,2,1)
plotSrvVals(2,2,1)

#--test in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate selectivity functions----
  info = inputs$lstSel;
  lstSelVals = calcSelectivity(dims,info,params,verbose);
  REPORT(lstSelVals);
  #--calculate survey catchability time series----
  info = inputs$lstSrv;
  lstSrvVals = calcFisheriesRates(dims,info,params,lstSelVals,verbose);
  REPORT(lstSrvVals);

  nll = -dnorm(1,params$dummy,1,log=TRUE);
  return(nll);
}

verbose=FALSE;
params$dummy = 0;
obj = MakeADFun(obj_fun,params,random=NULL,map=map,silent=FALSE);
rep = obj$report();
lstSrvVals = rep$lstSrv;
plotSrvVals(1,1,1)
plotSrvVals(2,1,1)
plotSrvVals(1,2,1)
plotSrvVals(2,2,1)

