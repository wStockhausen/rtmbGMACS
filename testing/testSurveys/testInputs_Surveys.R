#--test surveys specifications for the gmacs objective function
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
  source(file.path(dirPrj,"R","readParamInfo_Surveys.R"))
  source(file.path(dirPrj,"R","extractParamInfoFunctionType1.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Surveys.R"))
  source(file.path(dirPrj,"R","Functions_Catchbility.R"))
  source(file.path(dirPrj,"R","calcSurveys.R"))
}

##--set up model dimensions----
#----"y","s","r","x","m","a","z"
source(file.path(dirPrj,"testing/r_setupModelDimensions.TestA.R"))

type = "pre-specified-vertical";
if (type=="pre-specified-vertical"){
  ###--survey catchability with pre-specified-vertical format----
  dims = setupModelDims(zcs=seq(25.5,80.5,5));
  conn   = file.path(dirPrj,"testing/testSurveys/inputSpecs_Surveys.pre-specified-vertical.txt");
  res    = readParamInfo_Surveys(conn,TRUE);
  lstSrv = extractParamInfo_Surveys(res,dims,FALSE);
  params = list(pSrv_FPs=lstSrv$params);#--"FP" for "fixed" parameters
  map = lstSrv$map;
} else
if (type=="pre-specified-horizontal"){
  ###--survey catchability inputs with pre-specified-horizontal format----
  dims   = setupModelDims(zcs=seq(25.5,80.5,5));
  conn   = file.path(dirPrj,"testing/testSurveys/inputSpecs_Surveys.pre-specified-horizontal.txt");
  res    = readParamInfo_Surveys(conn,TRUE);
  lstSrv = extractParamInfo_Surveys(res,dims);
  params = list(pSrv_FPs=lstSrv$params);#--"FP" for "fixed" parameters
  map = lstSrv$map;
} else
if (type=="function"){
  dims = setupModelDims(zcs=seq(25.5,84.5,5));
  ###--selectivity functions with function type input----
  source(file.path(dirPrj,"R","readParamInfo_Selectivity.R"))
  source(file.path(dirPrj,"R","extractParamInfo_Selectivity.R"))
  source(file.path(dirPrj,"R","Functions_Selectivity.R"))
  source(file.path(dirPrj,"R","calcSelectivity.R"))
  conn   = file.path(dirPrj,"testing/testSelectivity/inputSpecs_Selectivity.function.txt");
  res    = readParamInfo_Selectivity(conn,FALSE);
  lstSel = extractParamInfo_Selectivity(res,dims,FALSE);
  params = list(pSel_MPs=lstSel$MPs$params);
  map = lstSel$map;
  if (!is.null(lstSel$OPs$params)) params[["pSel_OPs"]]=lstSel$OPs$params;
  if (!is.null(lstSel$DPs$params)) params[["pSel_DPs"]]=lstSel$DPs$params;
  if (!is.null(lstSel$REs$params)) params[["pSel_REs"]]=lstSel$REs$params;

  ###--survey catchability inputs with function format----
  conn   = file.path(dirPrj,"testing/testSurveys/inputSpecs_Surveys.function.txt");
  res    = readParamInfo_Surveys(conn,FALSE);
  lstSrv = extractParamInfo_Surveys(res,dims,FALSE);
  params = list(pSrv_MPs=lstSrv$MPs$params);
  map = lstSrv$map;
  if (!is.null(lstSrv$OPs$params)) params[["pSel_OPs"]]=lstSrv$OPs$params;
  if (!is.null(lstSrv$DPs$params)) params[["pSel_DPs"]]=lstSrv$DPs$params;
  if (!is.null(lstSrv$REs$params)) params[["pSel_REs"]]=lstSrv$REs$params;
}
inputs        = list();
inputs$dims   = dims;
inputs$lstSrv = lstSrv;#--add lstSrv to inputs

#--test calcSurveys function----
source(file.path(dirPrj,"R/calcSurveys.R"));
lstSrv = calcSurveys(inputs$dims,inputs$lstSrv,params,TRUE);
plotSrvVals<-function(if_=1,iy_=1,is_=1){
  dfrSel = dims$dmsC |>
            dplyr::mutate(value=(lstSrvVals[[if_]])[iy_,is_,],
                          z=as.numeric(as.character(z)));
  ggplot(dfrSrv, aes(x=as.numeric(z),y=value)) +
    geom_line() + geom_hline(yintercept=0,linetype=3,colour="white") +
    facet_grid(m+p~x) +
    labs(x="size",y="survey catchability") +
    wtsPlots::getStdTheme();
}
plotSrvVals(1,1,1)
plotSrvVals(2,1,1)
plotSrvVals(1,4,1)
plotSrvVals(2,4,1)

#--test arrSel in RTMB objective function----
obj_fun<-function(params){
  #--get dimensions----
  dims = inputs$dims;
  #--calculate fisheries time series----
  info = inputs$lstSrv;
  lstSrv = calcSurveys(dims,info,params,input$lstSel,verbose);
  REPORT(lstSrv);

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
plotSrvVals(1,4,1)
plotSrvVals(2,4,1)
